#include "custom_dynamics_world.h"

void CustomDynamicsWorld::internalSingleStepSimulation(btScalar timeStep)
{
    // TODO: Enable this? Don't think it's so important for our purposes here though.
    // BT_PROFILE("internalSingleStepSimulation");

    if (0 != m_internalPreTickCallback)
    {
        (*m_internalPreTickCallback)(this, timeStep);
    }

    btDispatcherInfo& dispatchInfo = getDispatchInfo();

    dispatchInfo.m_timeStep = timeStep;
    dispatchInfo.m_stepCount = 0;
    dispatchInfo.m_debugDraw = getDebugDrawer();

    getSolverInfo().m_timeStep = timeStep;

    integrateConstrainedBodiesWithCustomPhysics(timeStep);

    ///update vehicle simulation
    updateActions(timeStep);

    updateActivationState(timeStep);

    if (0 != m_internalTickCallback)
    {
        (*m_internalTickCallback)(this, timeStep);
    }
}

vector<btRigidBody *> collectBodies(btCollisionObjectArray & objects,
                                          vector<btRigidBody *> && buffer = vector<btRigidBody *>())
{
    const auto num_obj = objects.size();
    for (int i = 0; i < num_obj; ++i)
    {
        const auto & obj = objects[i];
        const auto body = btRigidBody::upcast(obj);

        if (body)
        {
            buffer.push_back(body);
        }
    }

    return buffer;
}

vector<btPersistentManifold*> fetchManifolds (btDynamicsWorld& world){
    const auto num_manifolds = world.getDispatcher()->getNumManifolds();
    vector<btPersistentManifold*> manifolds;
    for (int i = 0; i < num_manifolds; ++i) {
        manifolds.push_back(world.getDispatcher()->getManifoldByIndexInternal(i));
    }
    return manifolds;
}

/* sequential impulses method:
    1. Reset applied impulses α1 (and α2) and compute tangent vectors t1 (and t2) for each contact
    2. Until convergence or maximum number of iterations:
        2.1 Apply correctional impulses for joints
        2.2 Apply correctional normal impulses
        2.3 Apply correctional friction impulses */
void CustomDynamicsWorld::sequentialImpulses(btScalar timeStep){

    // filter out the btPoint2PointConstraint instances:
    const auto num_constraints = this->getNumConstraints();
    vector<btPoint2PointConstraint *> point_constraints;
    vector<btHingeConstraint *> hinge_constraints;
    vector<btScalar> accumulated_impulses;

    for (int i = 0; i < num_constraints; ++i){
        auto c = this->getConstraint(i);
        if(c->getConstraintType() == POINT2POINT_CONSTRAINT_TYPE){
            point_constraints.push_back(dynamic_cast<btPoint2PointConstraint *>(c));
        }else if(c->getConstraintType() == HINGE_CONSTRAINT_TYPE){
            hinge_constraints.push_back(dynamic_cast<btHingeConstraint *>(c));
            accumulated_impulses.push_back(0);
        }
    }

    //fetch all contact Manifolds
    auto manifolds = fetchManifolds(*this);
    vector<btScalar> target_velocities;

    //2
    for (int i = 0; i < getConstraintIterations(); i++){  
        //2.1 Ball joints constraints
        if(getApplyBallJointsCorrections()) point2PointConstraintCorrection(point_constraints, timeStep);

        if(getApplyHingeJointsCorrections()) hingeJointConstraintCorrection(hinge_constraints, timeStep, accumulated_impulses);
        // 2.2 Non-penetration and friction constraints
        if(getApplyFrictionCorrections() || getApplyContactCorrections()) manifoldCorrection(manifolds, timeStep, i, target_velocities);
    }
}

/*2.1 Foreach Constraint i: Compute Si for the system of the two rigid bodies constrained only by the current constraint in isolation
2.2 Compute an impulse represented by ∆λ˜_i by ∆λ˜_i = S^{−1]_i(−Giu), where u = (u_j , u_k).
2.3 Apply the impulse by updating the velocities of rigid bodies j and k, i.e. u <- u^p + M^{-1}G^T_i ∆λ˜*/
void CustomDynamicsWorld::point2PointConstraintCorrection(vector<btPoint2PointConstraint *> &constraints, btScalar timeStep){
    for (auto c : constraints){

        btRigidBody& body_j = c->getRigidBodyA();
        btRigidBody& body_k = c->getRigidBodyB();

        // 2.1 
        const auto r_j = c->getPivotInA(); // attachement point for body j
        const auto r_k = c->getPivotInB(); // attachement point for body k
        const auto R_j = body_j.getWorldTransform().getBasis();
        const auto R_k = body_k.getWorldTransform().getBasis();
        
        const auto m_inv_j = body_j.getInvMass();
        const auto m_inv_k = body_k.getInvMass();

        auto tensor_j = body_j.getInvInertiaTensorWorld();
        auto tensor_k = body_k.getInvInertiaTensorWorld();

        btVector3 v1,v2,v3;
        (R_j * r_j).getSkewSymmetricMatrix(&v1,&v2,&v3);
        btMatrix3x3 K_j = btMatrix3x3(v1,v2,v3);
        (R_k * r_k).getSkewSymmetricMatrix(&v1,&v2,&v3);
        btMatrix3x3 K_k = btMatrix3x3(v1,v2,v3);

        //Calculate S = GM^-1G^T
        btMatrix3x3 S = I*(m_inv_j+m_inv_k) - K_j*tensor_j*K_j - K_k*tensor_k*K_k;

        //2.2

        //Get Velocities
        btVector3 lV_j = body_j.getLinearVelocity();
        btVector3 aV_j = body_j.getAngularVelocity();
        btVector3 lV_k = body_k.getLinearVelocity();
        btVector3 aV_k = body_k.getAngularVelocity();

        //Calculate impulse
        btVector3 dC = lV_j-K_j*aV_j-lV_k+K_k*aV_k; // G*u
        
        btVector3 C = (body_j.getCenterOfMassPosition() + R_j*r_j)-(body_k.getCenterOfMassPosition() + R_k*r_k);

        btVector3 target_velocity = (-getGamma())*C*(1/timeStep) - dC;

        btVector3 impulse = {0,0,0};
        if(!isZero(S.determinant())){
            impulse = S.inverse() * target_velocity;
        }

        if(isZero(impulse)) continue;

        // 2.3
        body_j.setLinearVelocity(lV_j + m_inv_j*impulse);
        body_j.setAngularVelocity(aV_j + tensor_j*K_j*impulse);

        body_k.setLinearVelocity(lV_k - m_inv_k*impulse);
        body_k.setAngularVelocity(aV_k - tensor_k*K_k*impulse);
    }
}

/**/
void CustomDynamicsWorld::hingeJointConstraintCorrection(vector<btHingeConstraint *> &constraints, btScalar timeStep, vector<btScalar> &accumulated_impulses){
    for (int c_ind = 0; c_ind < constraints.size(); c_ind++){
        auto c = constraints.at(c_ind);

        hingeBallJointConstraint(c, timeStep);

        if(getHingeWith2x2()){
            hingeCombinedAxisConstraint(c, timeStep);
        }else{
            for(int i = 0; i < getHingeIterations(); i++){
                hingeIndividualAxisConstraint(c, timeStep);
            }
        }

        if(c->getEnableAngularMotor()) hingeMotorConstraint(c, c_ind, timeStep, accumulated_impulses);
    }
}

void CustomDynamicsWorld::hingeBallJointConstraint(btHingeConstraint* c, btScalar timeStep){
    btRigidBody& body_j = c->getRigidBodyA();
    btRigidBody& body_k = c->getRigidBodyB();

    const auto r_j = c->getFrameOffsetA().getOrigin(); // Ball joint connection point for body j
    const auto r_k = c->getFrameOffsetB().getOrigin(); // Ball joint connection point for body k
    const auto R_j = body_j.getWorldTransform().getBasis();
    const auto R_k = body_k.getWorldTransform().getBasis();
    
    const auto m_inv_j = body_j.getInvMass();
    const auto m_inv_k = body_k.getInvMass();

    auto tensor_j = body_j.getInvInertiaTensorWorld();
    auto tensor_k = body_k.getInvInertiaTensorWorld();

    btVector3 v1,v2,v3;
    (R_j * r_j).getSkewSymmetricMatrix(&v1,&v2,&v3);
    btMatrix3x3 K_j = btMatrix3x3(v1,v2,v3);
    (R_k * r_k).getSkewSymmetricMatrix(&v1,&v2,&v3);
    btMatrix3x3 K_k = btMatrix3x3(v1,v2,v3);

    //Calculate S = GM^-1G^T
    btMatrix3x3 S = I*(m_inv_j+m_inv_k) - K_j*tensor_j*K_j - K_k*tensor_k*K_k;

    //Get Velocities
    btVector3 lV_j = body_j.getLinearVelocity();
    btVector3 aV_j = body_j.getAngularVelocity();
    btVector3 lV_k = body_k.getLinearVelocity();
    btVector3 aV_k = body_k.getAngularVelocity();

    //Calculate impulse
    btVector3 dC = lV_j-K_j*aV_j-lV_k+K_k*aV_k; // G*u
    
    btVector3 C = (body_j.getCenterOfMassPosition() + R_j*r_j)-(body_k.getCenterOfMassPosition() + R_k*r_k);

    btVector3 target_velocity = (-getGamma())*C*(1/timeStep) - dC;

    btVector3 impulse = {0,0,0};
    if(!isZero(S.determinant())){
        impulse = S.inverse() * target_velocity;
    }

    if(!isZero(impulse)){
        body_j.setLinearVelocity(lV_j + m_inv_j*impulse);
        body_j.setAngularVelocity(aV_j + tensor_j*K_j*impulse);

        body_k.setLinearVelocity(lV_k - m_inv_k*impulse);
        body_k.setAngularVelocity(aV_k - tensor_k*K_k*impulse);
    }
}

void CustomDynamicsWorld::hingeCombinedAxisConstraint(btHingeConstraint* c, btScalar timeStep){
    btRigidBody& body_j = c->getRigidBodyA();
    btRigidBody& body_k = c->getRigidBodyB();

    const auto R_j = body_j.getWorldTransform().getBasis();
    const auto R_k = body_k.getWorldTransform().getBasis();

    // Bullet stores the orthogonal basis p, q, h as columns of an "offset frame" basis
    const auto h_j = c->getFrameOffsetA().getBasis().getColumn(2);
    const auto p_k = c->getFrameOffsetB().getBasis().getColumn(0);
    const auto q_k = c->getFrameOffsetB().getBasis().getColumn(1);

    auto tensor_j = body_j.getInvInertiaTensorWorld();
    auto tensor_k = body_k.getInvInertiaTensorWorld();

    btVector3 v1,v2,v3;
    (R_k * p_k).getSkewSymmetricMatrix(&v1,&v2,&v3);
    btMatrix3x3 K_p = btMatrix3x3(v1,v2,v3);
    (R_k * q_k).getSkewSymmetricMatrix(&v1,&v2,&v3);
    btMatrix3x3 K_q = btMatrix3x3(v1,v2,v3);
    btVector3 h_world = R_j * h_j;
    
    //Get Velocities
    btVector3 aV_j = body_j.getAngularVelocity();
    btVector3 aV_k = body_k.getAngularVelocity();

    cpMatrix dC(2,1);
    dC(0,0) = h_world.dot((K_p*(aV_j - aV_k)));
    dC(1,0) = h_world.dot((K_q*(aV_j - aV_k)));

    cpMatrix S(2,2);
    btMatrix3x3 tensor_sum = tensor_j+tensor_k;
    btMatrix3x3 temp = K_p*tensor_sum*K_p;
    S(0,0) = -multiplyVector3withMatrix3x3FromBothSides(h_world, temp);
    temp = K_p*tensor_sum*K_q;
    S(0,1) = -multiplyVector3withMatrix3x3FromBothSides(h_world, temp);
    temp = K_q*tensor_sum*K_p;
    S(1,0) = -multiplyVector3withMatrix3x3FromBothSides(h_world, temp);
    temp = K_q*tensor_sum*K_q;
    S(1,1) = -multiplyVector3withMatrix3x3FromBothSides(h_world, temp);

    cpMatrix impulse(2,1);
    if(!isZero(S.det2x2())){
        cpMatrix C(2,1);
        C(0,0) = h_world.dot(R_k*p_k);
        C(1,0) = h_world.dot(R_k*q_k);
        cpMatrix target_velocity = C*(1/timeStep)*getGamma()*(-1);
        impulse = S.invert2x2() * (target_velocity - dC);
    }

    if(!(isZero(impulse(0,0)) && isZero(impulse(1,0)))){
        body_j.setAngularVelocity(aV_j - tensor_j*K_p*h_world*impulse(0,0)
                                    - tensor_j*K_q*h_world*impulse(1,0));
        body_k.setAngularVelocity(aV_k + tensor_k*K_p*h_world*impulse(0,0)
                                    + tensor_k*K_q*h_world*impulse(1,0));
    }
}

void CustomDynamicsWorld::hingeIndividualAxisConstraint(btHingeConstraint* c, btScalar timeStep){
    btRigidBody& body_j = c->getRigidBodyA();
    btRigidBody& body_k = c->getRigidBodyB();

    const auto R_j = body_j.getWorldTransform().getBasis();
    const auto R_k = body_k.getWorldTransform().getBasis();

    // Bullet stores the orthogonal basis p, q, h as columns of an "offset frame" basis
    const auto h_j = c->getFrameOffsetA().getBasis().getColumn(2);
    const auto p_k = c->getFrameOffsetB().getBasis().getColumn(0);
    const auto q_k = c->getFrameOffsetB().getBasis().getColumn(1);

    auto tensor_j = body_j.getInvInertiaTensorWorld();
    auto tensor_k = body_k.getInvInertiaTensorWorld();

    btVector3 v1,v2,v3;
    (R_k * p_k).getSkewSymmetricMatrix(&v1,&v2,&v3);
    btMatrix3x3 K_p = btMatrix3x3(v1,v2,v3);
    (R_k * q_k).getSkewSymmetricMatrix(&v1,&v2,&v3);
    btMatrix3x3 K_q = btMatrix3x3(v1,v2,v3);
    btVector3 h_world = R_j * h_j;
    
    //Get Velocities
    btVector3 aV_j = body_j.getAngularVelocity();
    btVector3 aV_k = body_k.getAngularVelocity();

    btScalar dC_p = h_world.dot((K_p*(aV_j - aV_k)));
    btScalar dC_q = h_world.dot((K_q*(aV_j - aV_k)));

    btMatrix3x3 temp = K_p * (tensor_j + tensor_k) * K_p;
    btScalar S_p = -multiplyVector3withMatrix3x3FromBothSides(h_world, temp);
    temp = K_q * (tensor_j + tensor_k) * K_q;
    btScalar S_q = -multiplyVector3withMatrix3x3FromBothSides(h_world, temp);

    btScalar impulse_p = 0;
    if(!isZero(S_p)){
        btScalar C = h_world.dot(R_k*p_k);
        btScalar target_velocity = -getGamma()*C*(1/timeStep);
        impulse_p = (target_velocity-dC_p) * (1/S_p);
    }
    
    btScalar impulse_q = 0;
    if(!isZero(S_q)){
        btScalar C = h_world.dot(R_k*q_k);
        btScalar target_velocity = -getGamma()*C*(1/timeStep);
        impulse_q = (target_velocity-dC_q) * (1/S_q);
    }

    if(!(isZero(impulse_p) && isZero(impulse_q))){
        body_j.setAngularVelocity(aV_j - tensor_j*K_p*h_world*impulse_p
                                    - tensor_j*K_q*h_world*impulse_q);
        body_k.setAngularVelocity(aV_k + tensor_k*K_p*h_world*impulse_p
                                    + tensor_k*K_q*h_world*impulse_q);
    }
}

void CustomDynamicsWorld::hingeMotorConstraint(btHingeConstraint* c, int c_ind, btScalar timeStep, vector<btScalar> &accumulated_impulses){
    btRigidBody& body_j = c->getRigidBodyA();
    btRigidBody& body_k = c->getRigidBodyB();

    const auto R_j = body_j.getWorldTransform().getBasis();
    const auto R_k = body_k.getWorldTransform().getBasis();

    // Bullet stores the orthogonal basis p, q, h as columns of an "offset frame" basis
    const auto h_j = c->getFrameOffsetA().getBasis().getColumn(2);
    const auto h_k = c->getFrameOffsetB().getBasis().getColumn(2);

    auto tensor_j = body_j.getInvInertiaTensorWorld();
    auto tensor_k = body_k.getInvInertiaTensorWorld();

    btVector3 h = (R_j*h_j + R_k*h_k)/2;
    if(h.length() != 0) h.normalize();

    btVector3 aV_j = body_j.getAngularVelocity();
    btVector3 aV_k = body_k.getAngularVelocity();
    btScalar aV_target = c->getMotorTargetVelocity();

    btScalar dC = h.dot(aV_j-aV_k) - aV_target;

    btMatrix3x3 tensor_sum = tensor_j+tensor_k;
    btScalar S = multiplyVector3withMatrix3x3FromBothSides(h, tensor_sum);

    btScalar impulse = 0;
    if(!isZero(S)){
        impulse = (1/S) * -dC;

        auto& acc_impulse = accumulated_impulses.at(c_ind);
        auto max_impulse = c->getMaxMotorImpulse();

        btClamp<btScalar>(impulse, -1 * max_impulse - acc_impulse, max_impulse - acc_impulse);

        acc_impulse += impulse;
    }

    body_j.setAngularVelocity(aV_j + tensor_j*h*impulse);
    body_k.setAngularVelocity(aV_k - tensor_k*h*impulse);
}


/*
Combines contactCorrection and frictionCorrection as well as computation of t1, t2 to need to loop over all manifolds only once.
*/
void CustomDynamicsWorld::manifoldCorrection(vector<btPersistentManifold *> &manifolds, btScalar timeStep, int iteration, vector<btScalar> &target_velocities){
    auto num_manifolds = this->getDispatcher()->getNumManifolds();

    int tV_index = 0;
    for(size_t i = 0; i < num_manifolds; ++i){
        auto manifold = manifolds[i];

        // The const_cast is not nice, but apparently necessary
        auto body_j = btRigidBody::upcast(const_cast<btCollisionObject *>(manifold->getBody0()));
        auto body_k = btRigidBody::upcast(const_cast<btCollisionObject *>(manifold->getBody1()));

        const auto is_rigidJ = body_j != nullptr;
        const auto is_rigidK = body_k != nullptr;

        if(!is_rigidJ || !is_rigidK) {
            //Only do contact handling when both bodies are Rigidbodies
            continue;
        }

        const auto num_contacts = manifold->getNumContacts();
        for(int c = 0; c < num_contacts; ++c){
            auto& contact = manifold->getContactPoint(c);

            //Get all Variables
            auto n = contact.m_normalWorldOnB;

            const auto r_j = contact.m_localPointA;
            const auto r_k = contact.m_localPointB;
            auto R_j = body_j->getWorldTransform().getBasis();
            auto R_k = body_k->getWorldTransform().getBasis();

            const auto m_inv_j = body_j->getInvMass();
            const auto m_inv_k = body_k->getInvMass();

            auto tensor_j = body_j->getInvInertiaTensorWorld();
            auto tensor_k = body_k->getInvInertiaTensorWorld();

            btVector3 v1,v2,v3;
            (R_j * r_j).getSkewSymmetricMatrix(&v1,&v2,&v3);
            btMatrix3x3 K_j = btMatrix3x3(v1,v2,v3);
            (R_k * r_k).getSkewSymmetricMatrix(&v1,&v2,&v3);
            btMatrix3x3 K_k = btMatrix3x3(v1,v2,v3);

            //Get Velocities
            btVector3 lV_j = body_j->getLinearVelocity();
            btVector3 aV_j = body_j->getAngularVelocity();
            btVector3 lV_k = body_k->getLinearVelocity();
            btVector3 aV_k = body_k->getAngularVelocity();

            // if first iteration reset applied impulse α1 (and α2) and compute tangent vectors t1 (and t2) for each contact
            if (iteration == 0){                  
                // relative velocity at the contact point v_r = derivative C_b = G_b*u
                btVector3 v_r = lV_j-K_j*aV_j-lV_k+K_k*aV_k;

                // tangent velocity v_T is the projection of v_r onto the tangent plane
                btVector3 v_t = v_r - (v_r.dot(n))*n;

                if(v_t.length() == 0) contact.m_lateralFrictionDir1 = v_t;
                else contact.m_lateralFrictionDir1 = v_t / v_t.length();
                
                contact.m_lateralFrictionDir2 = n.cross(contact.m_lateralFrictionDir1);
                if(contact.m_lateralFrictionDir2.length() != 0) contact.m_lateralFrictionDir2.normalize();
                
                if (!getWarmStarting()){
                    contact.m_appliedImpulse = 0;
                    contact.m_appliedImpulseLateral1 = 0;
                    contact.m_appliedImpulseLateral2 = 0;
                }
            }
            
            auto KJK_j = K_j*tensor_j*K_j;
            auto KJK_k = K_k*tensor_k*K_k;

            // CONTACT CORRECTION:
            if(getApplyContactCorrections()){     
                btVector3 worldDiff = (body_j->getCenterOfMassPosition() + R_j*r_j)-(body_k->getCenterOfMassPosition() + R_k*r_k);
                btScalar C =  n.dot(worldDiff);
                btScalar dC = n.dot(lV_j)-n.dot(K_j*aV_j)-n.dot(lV_k)+n.dot(K_k*aV_k); //G*u 

                if(iteration == 0){
                    btScalar stabilization = (-getGamma())*C*(1/timeStep);
                    btScalar restitution = (-contact.m_combinedRestitution) * dC;
                    btScalar tV = max(stabilization, restitution);
                    if(C > epsilon) tV = 0;
                    target_velocities.push_back(tV);
                    //cout << contact.m_targetVelocity << endl;
                }

                if(C < -epsilon || (isZero(C) && dC < -epsilon)){ //Constraint is only violated when C negative or when C = 0 but dC<0

                    btScalar impulse = 0;
                    if(iteration == 0 && getWarmStarting()){
                        contact.m_appliedImpulse *= getWarmStartingFactor();
                        impulse = contact.m_appliedImpulse;
                    }else{
                        //We only need this S when we reach this conditional
                        btScalar S = (m_inv_j + m_inv_k)*n.length2() - multiplyVector3withMatrix3x3FromBothSides(n, KJK_j)
                                                                    - multiplyVector3withMatrix3x3FromBothSides(n, KJK_k);
                        if(!isZero(S)){ // S needs to be invertible
                            impulse = (target_velocities.at(tV_index+c) - dC) * (1/S);
                        }

                        if(isZero(impulse)) continue;

                        if(contact.m_appliedImpulse + impulse < 0) {
                            impulse = -contact.m_appliedImpulse;
                        }
                        contact.m_appliedImpulse += impulse;
                    }
    
                    //apply impulse
                    body_j->setLinearVelocity(lV_j+m_inv_j*impulse*n);
                    body_j->setAngularVelocity(aV_j+impulse*(tensor_j*K_j*n));

                    body_k->setLinearVelocity(lV_k-m_inv_k*impulse*n);
                    body_k->setAngularVelocity(aV_k-impulse*(tensor_k*K_k*n));
                }
            }

            // FRICTION CORRECTION:
            if(getApplyFrictionCorrections()) {
                auto t1 = contact.m_lateralFrictionDir1;
                auto t2 = contact.m_lateralFrictionDir2;

                R_j = body_j->getWorldTransform().getBasis();
                R_k = body_k->getWorldTransform().getBasis();

                tensor_j = body_j->getInvInertiaTensorWorld();
                tensor_k = body_k->getInvInertiaTensorWorld();

                (R_j * r_j).getSkewSymmetricMatrix(&v1,&v2,&v3);
                K_j = btMatrix3x3(v1,v2,v3);
                (R_k * r_k).getSkewSymmetricMatrix(&v1,&v2,&v3);
                K_k = btMatrix3x3(v1,v2,v3);

                KJK_j = K_j*tensor_j*K_j;
                KJK_k = K_k*tensor_k*K_k;
                
                // get new velocities, since we may have changed them for contact correction
                lV_j = body_j->getLinearVelocity();
                aV_j = body_j->getAngularVelocity();
                lV_k = body_k->getLinearVelocity();
                aV_k = body_k->getAngularVelocity();

                auto mu = getMU();
                auto lambda = contact.m_appliedImpulse;

                btScalar deltaA1 = 0;
                btScalar deltaA2 = 0;
                if(iteration == 0 && getWarmStarting()){
                    contact.m_appliedImpulseLateral1 *= getWarmStartingFactor();
                    contact.m_appliedImpulseLateral2 *= getWarmStartingFactor();
                    deltaA1 = contact.m_appliedImpulseLateral1;
                    deltaA2 = contact.m_appliedImpulseLateral2;
                }else{
                    //Calculate S = GM^-1G^T
                    btScalar S1 = (m_inv_j + m_inv_k)*t1.length2() - multiplyVector3withMatrix3x3FromBothSides(t1, KJK_j)
                                                            - multiplyVector3withMatrix3x3FromBothSides(t1, KJK_k);
                    btScalar S2 = (m_inv_j + m_inv_k)*t2.length2() - multiplyVector3withMatrix3x3FromBothSides(t2, KJK_j)
                                                            - multiplyVector3withMatrix3x3FromBothSides(t2, KJK_k);

                    if(!isZero(S1)){
                        auto& a1 = contact.m_appliedImpulseLateral1;

                        btScalar dC1 = t1.dot(lV_j)-t1.dot(K_j*aV_j)-t1.dot(lV_k)+t1.dot(K_k*aV_k);
                        deltaA1 = dC1 * -1 * (1/S1);

                        btClamp<btScalar>(deltaA1, -1 * (mu * lambda) - a1, (mu * lambda) - a1);

                        if(isZero(deltaA1)) deltaA1 = 0;
                        a1 += deltaA1;
                    }
                    if(!isZero(S2)){        
                        auto& a2 = contact.m_appliedImpulseLateral2;

                        btScalar dC2 = t2.dot(lV_j)-t2.dot(K_j*aV_j)-t2.dot(lV_k)+t2.dot(K_k*aV_k);
                        deltaA2 = dC2 * -1 * (1/S2);

                        btClamp<btScalar>(deltaA2, -1 * (mu * lambda) - a2, (mu * lambda) - a2);

                        if(isZero(deltaA2)) deltaA2 = 0;
                        a2 += deltaA2;
                    }
                }

                //apply impulses
                body_j->setLinearVelocity(lV_j+m_inv_j*(deltaA1*t1 + deltaA2*t2));
                body_j->setAngularVelocity(aV_j+tensor_j*(K_j*t1*deltaA1 + K_j*t2*deltaA2));

                body_k->setLinearVelocity(lV_k-m_inv_k*(deltaA1*t1 + deltaA2*t2));
                body_k->setAngularVelocity(aV_k-tensor_k*(K_k*t1*deltaA1 + K_k*t2*deltaA2));
            }
        }

        tV_index += manifold->getNumContacts();
    }
}

void CustomDynamicsWorld::integrateConstrainedBodiesWithCustomPhysics(btScalar timeStep) {
    // We typically perform collision detection at the beginning of the step
    performDiscreteCollisionDetection();

    // Collect all instances of btRigidBody (as pointers) in the current simulation. Note: there may also be
    // collision objects which are not rigid bodies.
    auto bodies = collectBodies(getCollisionObjectArray());

    //Update velocity for all bodies (for constrained bodies we calculate u_p here)
    for(auto body : bodies){
        body->setLinearVelocity(body->getLinearVelocity() + timeStep * body->getInvMass() * body->getTotalForce());
        body->setAngularVelocity(body->getAngularVelocity() + body->getInvInertiaTensorWorld() * timeStep * body->getTotalTorque());
    }

    sequentialImpulses(timeStep);

    for (auto body : bodies)
    {
        auto x = body->getCenterOfMassPosition();
        auto q = body->getOrientation();
        body->applyGravity();
        
        x += timeStep * body->getLinearVelocity();
        btVector3 omega = body->getAngularVelocity();
        q += btQuaternion(omega.x(), omega.y(), omega.z(), 0) * q * timeStep * btScalar(0.5f);
        q.safeNormalize();
        body->setCenterOfMassTransform(btTransform(q, x));
    }

    // When looping over your constraints, you can use getConstraintIterations() to obtain the number of
    // constraint iterations configured in the Project Settings:
    //  for (int i = 0; i < getConstraintIterations(); ++i)
    //  {
    //    // process constraints
    //  }

    // Important types that you may wish to use:
    //  btVector3 (3-element vector)
    //  btMatrix3 (3x3 matrix)
    //  btTransform (encapsulates rotation transform and translation)

    // Matrix/vector types support "obvious" functionality, such as:
    //  v.cross(u)
    //  v.dot(u)
    //  a * v (matrix-vector product)
    //  a * b (matrix-matrix product)
    //  a * alpha (scaling matrix by scalar)

    // Below we have included some examples of things you may want to do
    // Given a rigid body `body`:
    //  const auto x = body->getCenterOfMassPosition();
    //  const auto q = body->getOrientation();
    //  const auto v = body->getLinearVelocity();
    //  const auto omega = body->getAngularVelocity();
    //  const auto m_inv = body->getInvMass();
    //  const auto I_inv = body->getInvInertiaTensorWorld();
    //  const auto f = body->getTotalForce();
    //  const auto tau = body->getTotalTorque();

    // Create quaternion from angular velocity (with 0 "w" component)
    //  const auto omega_quat = btQuaternion(omega.x(), omega.y(), omega.z(), 0);

    // Normalize a quaternion:
    //  q.safeNormalize();

    // Update transform (generalized position) of rigid body with new values q and x
    // (note: implicitly also updates "world" inertia tensor)
    //  body->setCenterOfMassTransform(btTransform(q, x));

    // Get rotation matrix for a body
    //  body->getWorldTransform().getBasis()
}

// =========================== HELPER FUNCTIONS ===========================

#pragma region Helper_Functions
/* calculates the constraint velocity map for a ball joint connecting two rigid bodies.

Parameters
R_j         btMatrix3x3 object representing the rotation matrix of the first body.
R_k		    btMatrix3x3 object representing the rotation matrix of the first body.
r_j         btVector3 object representing the attachment point of the first body.
r_k         btVector3 object representing the attachment point of the second body.

Returns
cpMatrix object representing the constraint velocity map G. It is a 3x12 matrix that 
describes the relationship between the velocities of the two rigid bodies connected by a ball joint.*/
cpMatrix CustomDynamicsWorld::computeConstraintVelocityMap(btMatrix3x3 R_j, btMatrix3x3 R_k,btVector3 r_j, btVector3 r_k){
    btVector3 v1,v2,v3;
    (R_j * r_j).getSkewSymmetricMatrix(&v1,&v2,&v3);
    btMatrix3x3 K_j = btMatrix3x3(v1,v2,v3);
    (R_k * r_k).getSkewSymmetricMatrix(&v1,&v2,&v3);
    btMatrix3x3 K_k = btMatrix3x3(v1,v2,v3);

    cpMatrix G(3,12);
    G.initializeBlock(I, 0, 0);
    G.initializeBlock(K_j*-1, 0, 3);
    G.initializeBlock(I*-1, 0, 6);
    G.initializeBlock(K_k, 0, 9);

    return G;
}

/* calculates the inverse mass matrix for a pair of rigid bodies.

Parameters
body_j      btRigidBody object representing the first rigid body.
body_k      btRigidBody object representing the second rigid body.

Returns
cpMatrix object representing the inverse mass matrix M_inv. It is a 12x12 matrix that 
is comprised by the inverse mass and inverse inertia tensors of the two rigid bodies.*/
cpMatrix CustomDynamicsWorld::computeInverseMassMatrix(btRigidBody body_j, btRigidBody body_k){
    btMatrix3x3 mI_j = I * body_j.getInvMass();
    btMatrix3x3 mI_k = I * body_k.getInvMass();
    btMatrix3x3 tensor_j = body_j.getInvInertiaTensorWorld();
    btMatrix3x3 tensor_k = body_k.getInvInertiaTensorWorld();

    // the block-diagonal matrix M = diag(Mj,Mk) // see crash_course.pdf 2.2 end, page 14
    cpMatrix M_inv(12,12); 
    M_inv.initializeBlock(mI_j, 0, 0);
    M_inv.initializeBlock(tensor_j, 3, 3);
    M_inv.initializeBlock(mI_k, 6, 6);
    M_inv.initializeBlock(tensor_k, 9, 9);

    return M_inv;
}


/* calculates a 12x1 matrix that represents the linear and angular velocities of the bodies.

Parameters
body_j      btRigidBody object representing the first rigid body.
body_k      btRigidBody object representing the second rigid body.

Returns
cpMatrix object representing the constraint velocity matrix u_mat. It is a 12x1 matrix that 
is comprised by the linear and angular velocities of the two rigid bodies.*/
cpMatrix CustomDynamicsWorld::computeConstraintVelocityMatrix(btRigidBody body_j, btRigidBody body_k){
    // Retrieve the linear and angular velocities of the rigid bodies
    btVector3 velocities[] = {body_j.getLinearVelocity(), body_j.getAngularVelocity(), body_k.getLinearVelocity(), body_k.getAngularVelocity()};
    // construct velocity matrix
    cpMatrix u_mat(12,1); //12x1
    for(int i = 0; i < 4; i++){
        u_mat.setWithBtVector3(velocities[i], i*3, 0);
    }

    return u_mat;
}

/* modifies the linear and angular velocities of the bodies based on the given impulse vector.

Parameters
body_j      btRigidBody object representing the first rigid body.
body_k      btRigidBody object representing the second rigid body.
impulse     cpMatrix object representing the impulse vector to be applied. The dimensions of impulse should be 12x1.*/
void CustomDynamicsWorld::applyImpulse(btRigidBody& body_j, btRigidBody& body_k, cpMatrix impulse){
    body_j.setLinearVelocity( body_j.getLinearVelocity()  + impulse.getBtVector3(0,0));
    body_j.setAngularVelocity(body_j.getAngularVelocity() + impulse.getBtVector3(3,0));
    body_k.setLinearVelocity( body_k.getLinearVelocity()  + impulse.getBtVector3(6,0));
    body_k.setAngularVelocity(body_k.getAngularVelocity() + impulse.getBtVector3(9,0));
}

/*calculates vec^T*mat*vec which isn't easily possible with only bullet datatypes

Parameters
vec      btVector3 that gets multiplied to the matrix.
mat      btMatrix3x3 its obvious isnt it.

Returns
btScalar, since the dimensions are 1x3 * 3x3 * 3*1 = 1*1.*/
btScalar CustomDynamicsWorld::multiplyVector3withMatrix3x3FromBothSides(btVector3& vec, btMatrix3x3& mat){
    btScalar res = mat[0][0]*vec.getX()*vec.getX() + mat[1][1]*vec.getY()*vec.getY() + mat[2][2]*vec.getZ()*vec.getZ()
                 + (mat[1][0]+mat[0][1])*vec.getX()*vec.getY()
                 + (mat[2][0]+mat[0][2])*vec.getX()*vec.getZ()
                 + (mat[2][1]+mat[1][2])*vec.getY()*vec.getZ();
    return res;
}

#pragma endregion Helper_Functions
#pragma region Prints

void CustomDynamicsWorld::printVector(const btVector3& vector, string name) {
    cout << "Vector " << name << ":"; 
    for (int i = 0; i < 3; i++) {
        cout << dec << vector[i] << ", ";
    }
    cout << endl;
}

void CustomDynamicsWorld::printMatrix(const btMatrix3x3& matrix, string name) {
    cout << "Matrix " << name << ":" << endl; 
    for (int row = 0; row < 3; row++) {
        for (int col = 0; col < 3; col++) {
            cout << dec << matrix[row][col] << ", ";
        }
        cout << endl;
    }
    cout << endl;
}

void CustomDynamicsWorld::printQuat(const btQuaternion quat, string name) {
    cout << name << ": ("; 
    cout << dec << quat.w() << ", ";
    cout << dec << quat.x() << ", ";
    cout << dec << quat.y() << ", ";
    cout << dec << quat.z() << ")\n";
}

#pragma endregion Prints