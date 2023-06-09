#include "custom_dynamics_world.h"

typedef btVector3 btVector12[4];

using std::cout;
using std::vector;
using std::string;
using std::endl;
using std::dec;

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

std::vector<btRigidBody *> collectBodies(btCollisionObjectArray & objects,
                                          std::vector<btRigidBody *> && buffer = std::vector<btRigidBody *>())
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

std::vector<btPersistentManifold*> fetchManifolds (btDynamicsWorld& world){
    const auto num_manifolds = world.getDispatcher()->getNumManifolds();
    std::vector<btPersistentManifold*> manifolds;
    for (int i = 0; i < num_manifolds; ++i) {
        manifolds.push_back(world.getDispatcher()->getManifoldByIndexInternal(i));
    }
    return manifolds;
}

#pragma region Prints

void printVector(const btVector3& vector, char name) {
    std::cout << name << ":"; 
  for (int i = 0; i < 3; i++) {
     std::cout << std::dec << vector[i] << ", ";
  }
  std::cout << std::endl;
}

void printMatrix(const btMatrix3x3& matrix, char name) {
    std::cout << name << ":"; 
    for (int row = 0; row < 3; row++) {
        for (int col = 0; col < 3; col++) {
            std::cout << std::dec << matrix[row][col] << ", ";
        }
        std::cout << std::endl;
    }
}

void printQuat(const btQuaternion quat, char name) {
    std::cout << name << ": ("; 
    std::cout << std::dec << quat.w() << ", ";
    std::cout << std::dec << quat.x() << ", ";
    std::cout << std::dec << quat.y() << ", ";
    std::cout << std::dec << quat.z() << ")\n";
}

#pragma endregion Prints

void computeOrthogonalVectors(const btVector3& n, btVector3& t1, btVector3& t2) {
    // Normalize n
    float length = std::sqrt(n.getX() * n.getX() + n.getY() * n.getY() + n.getZ() * n.getZ());
    btVector3 normalizedN = { n.getX() / length, n.getY() / length, n.getZ() / length };

    // Find a vector orthogonal to n
    if (normalizedN.getX() != 0 || normalizedN.getY() != 0) {
        t1 = { -normalizedN.getY(), normalizedN.getX(), 0 };
    } else {
        t1 = { 1, 0, 0 };
    }

    // Compute b by taking the cross product of n and a
    t2 = normalizedN.cross(t1);
}

/* sequential impulses method:
    2. Until convergence or maximum number of iterations:
        2.1 Apply correctional impulses for joints
        2.2 Apply correctional impulses for contacts */
void CustomDynamicsWorld::sequentialImpulses(btScalar timeStep){

    // filter out the btPoint2PointConstraint instances:
    const auto num_constraints = this->getNumConstraints();
    std::vector<btPoint2PointConstraint *> point_constraints;

    for (int i = 0; i < num_constraints; ++i){
        auto c = this->getConstraint(i);
        if(c->getConstraintType() == POINT2POINT_CONSTRAINT_TYPE){
            point_constraints.push_back(dynamic_cast<btPoint2PointConstraint *>(c));
        }
    }

    //fetch all contact Manifolds
    auto manifolds = fetchManifolds(*this);
    //reset applied_impulse
    for(size_t i = 0; i < this->getDispatcher()->getNumManifolds(); ++i){
        auto manifold = manifolds[i];
        
        auto body_j = btRigidBody::upcast((manifold->getBody0()));
        auto body_k = btRigidBody::upcast((manifold->getBody1()));

        const auto is_rigidJ = body_j != nullptr;
        const auto is_rigidK = body_k != nullptr;

        if(!is_rigidJ || !is_rigidK) {
            //Only do contact handling when both bodies are Rigidbodies
            continue;
        }

        const auto num_contacts = manifold->getNumContacts();
        for(int c = 0; c < num_contacts; ++c){
            auto& contact = manifold->getContactPoint(c);
            const auto n = contact.m_normalWorldOnB;

            const auto r_j = contact.m_localPointA;
            const auto r_k = contact.m_localPointB;

            contact.m_appliedImpulse = 0;
            contact.m_appliedImpulseLateral1 = 0;
            contact.m_appliedImpulseLateral2 = 0;

            const btMatrix3x3 I = btMatrix3x3::getIdentity();
            btMatrix3x3 R_j = body_j->getWorldTransform().getBasis();
            btMatrix3x3 R_k = body_k->getWorldTransform().getBasis();

            //printMatrix(R_j, 'R');

            //Calculate G
            btVector3 v1,v2,v3;
            (R_j * r_j).getSkewSymmetricMatrix(&v1,&v2,&v3);
            btMatrix3x3 K_j = btMatrix3x3(v1,v2,v3);
            (R_k * r_k).getSkewSymmetricMatrix(&v1,&v2,&v3);
            btMatrix3x3 K_k = btMatrix3x3(v1,v2,v3);

            cpMatrix G_b(3,12); //Ball Joint G
            G_b.initializeBlock(I, 0, 0);
            G_b.initializeBlock(K_j*-1, 0, 3);
            G_b.initializeBlock(I*-1, 0, 6);
            G_b.initializeBlock(K_k, 0, 9);

            //G_b.print("G_b");

            btVector12 u = {body_j->getLinearVelocity(), body_j->getAngularVelocity(), body_k->getLinearVelocity(), body_k->getAngularVelocity()};
            cpMatrix u_mat(12,1); //12x1
            for(int i = 0; i < 4; i++){
                u_mat.setWithBtVector3(u[i], i*3, 0);
            }

            //u_mat.print("u");

            // relative velocity at the contact point v_r = derivative C_b
            btVector3 v_r = (G_b*u_mat).toBtVector3();

            //printVector(v_r, 'C');

            // tangent velocity v_T is the projection of v_r onto the tangent plane
            btVector3 v_t = v_r - (v_r.dot(n))*n;
            
            //cout << "n . v_t: " << dec << n.dot(v_t) << endl;

            //assert(n.dot(v_t) == 0);

            //printVector(v_t, 'v');

            if(v_t.length() == 0) contact.m_lateralFrictionDir1 = {0,0,0};
            else contact.m_lateralFrictionDir1 = v_t / v_t.length();
            
            contact.m_lateralFrictionDir2 = n.cross(contact.m_lateralFrictionDir1);
            //printVector(contact.m_lateralFrictionDir2, 't');

            if(contact.m_lateralFrictionDir2.length() != 0) contact.m_lateralFrictionDir2.normalize();
            
            //printVector(contact.m_lateralFrictionDir2, 'a');

            
            //cout << "n . t2: " << dec << n.dot(v_t) << endl;
            //cout << "t1 . t2: " << dec << n.dot(v_t) << endl;

            //assert(n.dot(contact.m_lateralFrictionDir1) == 0);
            //assert(n.dot(contact.m_lateralFrictionDir2) == 0);
            //assert(contact.m_lateralFrictionDir1.dot(contact.m_lateralFrictionDir2) == 0);
        }
    }

    //2
    for (int i = 0; i < getConstraintIterations(); i++){  
        //2.1 Ball joint constraints
        point2PointConstraintCorrection(point_constraints, timeStep);

        //2.2 Non-penetration constraints
        contactCorrection(manifolds, timeStep);

        //2.3 Coulomb friction
        frictionCorrection(manifolds, timeStep);
    }
}

/*2.1 Foreach Constraint i: Compute Si for the system of the two rigid bodies constrained only by the current constraint in isolation
    2.1.1 Compute constraint velocity map G
    2.1.2 getting mass matrix M (or rather M^{-1})
    2.1.3 compute S
2.2 Compute an impulse represented by ∆λ˜_i by ∆λ˜_i = S^{−1]_i(−Giu), where u = (u_j , u_k).
2.3 Apply the impulse by updating the velocities of rigid bodies j and k, i.e. u <- u^p + M^{-1}G^T_i ∆λ˜
    (2.3.1 compute u^p = u^n + ∆t  M−1 f^{ext})
    2.3.2 compute M^{-1}G^T_i ∆λ
    2.3.3 compute new velocity u*/
void CustomDynamicsWorld::point2PointConstraintCorrection(std::vector<btPoint2PointConstraint *> &constraints, btScalar timeStep){
    for (auto c : constraints){

        btRigidBody& body_j = c->getRigidBodyA();
        btRigidBody& body_k = c->getRigidBodyB();

        const btVector3& r_j = c->getPivotInA(); // attachement point for body j
        const btVector3& r_k = c->getPivotInB(); // attachement point for body k
        const btMatrix3x3 I = btMatrix3x3::getIdentity();

        // 2.1 
        btMatrix3x3 R_j = body_j.getWorldTransform().getBasis();
        btMatrix3x3 R_k = body_k.getWorldTransform().getBasis();

        // 2.1.1
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
        cpMatrix G_transpose = G.transpose();

        // 2.1.2
        btMatrix3x3 e; // null matrix
        btMatrix3x3 mI_j = I * body_j.getInvMass();
        btMatrix3x3 mI_k = I * body_k.getInvMass();
        btMatrix3x3 tensor_j = body_j.getInvInertiaTensorWorld();
        btMatrix3x3 tensor_k = body_k.getInvInertiaTensorWorld();

        // the block-diagonal matrix M = diag(Mj,Mk) // see crash_course.pdf 2.2 end, page 14
        cpMatrix M_inv(12,12); // inverse mass matrix
        M_inv.initializeBlock(mI_j, 0, 0);
        M_inv.initializeBlock(tensor_j, 3, 3);
        M_inv.initializeBlock(mI_k, 6, 6);
        M_inv.initializeBlock(tensor_k, 9, 9);

        // 2.1.3
        btMatrix3x3 S = (G * M_inv * G_transpose).toBtMatrix3x3();

        btVector12 u = {body_j.getLinearVelocity(), body_j.getAngularVelocity(), body_k.getLinearVelocity(), body_k.getAngularVelocity()};
        cpMatrix u_mat(12,1);
        for(int i = 0; i < 4; i++){
            u_mat.setWithBtVector3(u[i], i*3, 0);
        }

        //2.2
        btVector3 ndC = (G * -1 * u_mat).toBtVector3();
        
        btVector3 C = (body_j.getCenterOfMassPosition() + R_j*r_j)-(body_k.getCenterOfMassPosition() + R_k*r_k);

        btVector3 target_veclotiy = (-getGamma())*C*(1/timeStep) + ndC;
        cpMatrix impulse(3,1);
        if(S.determinant() != 0){
            impulse.setWithBtVector3(S.inverse() * target_veclotiy, 0, 0);
        }

        // 2.3
        // 2.3.2
        cpMatrix M_invG_transposeDeltaLambda = M_inv * G_transpose * impulse; // 12x1
        // 2.3.3
        body_j.setLinearVelocity( u[0] + M_invG_transposeDeltaLambda.getBtVector3(0,0));
        body_j.setAngularVelocity(u[1] + M_invG_transposeDeltaLambda.getBtVector3(3,0));
        body_k.setLinearVelocity( u[2] + M_invG_transposeDeltaLambda.getBtVector3(6,0));
        body_k.setAngularVelocity(u[3] + M_invG_transposeDeltaLambda.getBtVector3(9,0));
    }
}

void CustomDynamicsWorld::contactCorrection(std::vector<btPersistentManifold *> &manifolds, btScalar timeStep){
    auto num_manifolds = this->getDispatcher()->getNumManifolds();

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
            const auto n = contact.m_normalWorldOnB;
            const auto r_j = contact.m_localPointA;
            const auto r_k = contact.m_localPointB;

            //Needed Variables
            cpMatrix n_cp(3,1);
            n_cp.setWithBtVector3(n, 0, 0);
            const btMatrix3x3 I = btMatrix3x3::getIdentity();
            btMatrix3x3 R_j = body_j->getWorldTransform().getBasis();
            btMatrix3x3 R_k = body_k->getWorldTransform().getBasis();

            //Compute and apply delta_lambda etc.

            //Calculate G
            btVector3 v1,v2,v3;
            (R_j * r_j).getSkewSymmetricMatrix(&v1,&v2,&v3);
            btMatrix3x3 K_j = btMatrix3x3(v1,v2,v3);
            (R_k * r_k).getSkewSymmetricMatrix(&v1,&v2,&v3);
            btMatrix3x3 K_k = btMatrix3x3(v1,v2,v3);

            cpMatrix G_b(3,12); //Ball Joint G
            G_b.initializeBlock(I, 0, 0);
            G_b.initializeBlock(K_j*-1, 0, 3);
            G_b.initializeBlock(I*-1, 0, 6);
            G_b.initializeBlock(K_k, 0, 9);

            cpMatrix G = n_cp.transpose() * G_b; // 1x12


            //Calculate S= GM^-1G^T
            btMatrix3x3 mI_j = I * body_j->getInvMass();
            btMatrix3x3 mI_k = I * body_k->getInvMass();
            btMatrix3x3 tensor_j = body_j->getInvInertiaTensorWorld();
            btMatrix3x3 tensor_k = body_k->getInvInertiaTensorWorld();

            // the block-diagonal matrix M = diag(Mj,Mk) // see crash_course.pdf 2.2 end, page 14
            cpMatrix M_inv(12,12); // inverse mass matrix
            M_inv.initializeBlock(mI_j, 0, 0);
            M_inv.initializeBlock(tensor_j, 3, 3);
            M_inv.initializeBlock(mI_k, 6, 6);
            M_inv.initializeBlock(tensor_k, 9, 9);

            btScalar S = (G*M_inv*G.transpose())(0,0); // The Product returns a 1x1 Matrix, which is why we just save S as a scalar

            //Calculate Impulse
            btVector12 u = {body_j->getLinearVelocity(), body_j->getAngularVelocity(), body_k->getLinearVelocity(), body_k->getAngularVelocity()};
            cpMatrix u_mat(12,1); //12x1
            for(int i = 0; i < 4; i++){
                u_mat.setWithBtVector3(u[i], i*3, 0);
            }

            cpMatrix worldDiff(3,1);
            worldDiff.setWithBtVector3(((body_j->getCenterOfMassPosition() + R_j*r_j)-(body_k->getCenterOfMassPosition() + R_k*r_k)), 0, 0);
            btScalar C =  (n_cp.transpose() * worldDiff)(0,0);
            btScalar dC = (G * u_mat)(0,0); //G * u_mat is 1x1 (1x12 * 12x1 = 1x1)

            if(C < 0 || (C==0 && dC<0)){ //Constraint is only violated when C negative or when C = 0 but dC<0
            
                // This is analogous to ball joints; To test this let a box fall onto a plane from high up 
                // => With Gamma = 0 it will clip into the plane but with Gamma > 0 it will correct itself
                // A nice effect is: The higher gamma the more the objects bounce, with Gamma > 1 they will gain energy with each collision
                btScalar target_veclotiy = (-getGamma())*C*(1/timeStep) - dC;

                btScalar impulse = 0;
                if(S != 0){ // S needs to be invertible and dont apply impulse if the bodies are moving apart (dc > 0)
                    impulse = target_veclotiy * (1/S);
                }

                if(contact.m_appliedImpulse + impulse < 0) continue;
                contact.m_appliedImpulse += impulse;

                //apply impulse
                cpMatrix M_invG_transposeDeltaLambda = M_inv * G.transpose() * impulse;

                body_j->setLinearVelocity( u[0] + M_invG_transposeDeltaLambda.getBtVector3(0,0));
                body_j->setAngularVelocity(u[1] + M_invG_transposeDeltaLambda.getBtVector3(3,0));
                body_k->setLinearVelocity( u[2] + M_invG_transposeDeltaLambda.getBtVector3(6,0));
                body_k->setAngularVelocity(u[3] + M_invG_transposeDeltaLambda.getBtVector3(9,0));
                
            }

        }
    }
}


void CustomDynamicsWorld::frictionCorrection(std::vector<btPersistentManifold *> &manifolds, btScalar timeStep){
    auto num_manifolds = this->getDispatcher()->getNumManifolds();

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
            const auto n = contact.m_normalWorldOnB;
            const auto t1 = contact.m_lateralFrictionDir1;
            const auto t2 = contact.m_lateralFrictionDir2;
            const auto r_j = contact.m_localPointA;
            const auto r_k = contact.m_localPointB;

            //Needed Variables
            cpMatrix n_cp(3,1);
            n_cp.setWithBtVector3(n, 0, 0);
            cpMatrix t1_cp(3,1);
            t1_cp.setWithBtVector3(t1, 0, 0);
            cpMatrix t2_cp(3,1);
            t2_cp.setWithBtVector3(t2, 0, 0);

            // t1_cp.print("t1");
            // t2_cp.print("t2");

            //cout << "a2: " << dec << a2 << endl;

            const btMatrix3x3 I = btMatrix3x3::getIdentity();
            btMatrix3x3 R_j = body_j->getWorldTransform().getBasis();
            btMatrix3x3 R_k = body_k->getWorldTransform().getBasis();

            //Compute and apply delta_lambda etc.

            //Calculate G
            btVector3 v1,v2,v3;
            (R_j * r_j).getSkewSymmetricMatrix(&v1,&v2,&v3);
            btMatrix3x3 K_j = btMatrix3x3(v1,v2,v3);
            (R_k * r_k).getSkewSymmetricMatrix(&v1,&v2,&v3);
            btMatrix3x3 K_k = btMatrix3x3(v1,v2,v3);

            cpMatrix G_b(3,12); //Ball Joint G
            G_b.initializeBlock(I, 0, 0);
            G_b.initializeBlock(K_j*-1, 0, 3);
            G_b.initializeBlock(I*-1, 0, 6);
            G_b.initializeBlock(K_k, 0, 9);

            //G_b.print("G_b");

            cpMatrix G1 = t1_cp.transpose() * G_b; // 1x12

            //G1.print("G1");

            cpMatrix G2 = t2_cp.transpose() * G_b;

            //Calculate S= GM^-1G^T
            btMatrix3x3 mI_j = I * body_j->getInvMass();
            btMatrix3x3 mI_k = I * body_k->getInvMass();
            btMatrix3x3 tensor_j = body_j->getInvInertiaTensorWorld();
            btMatrix3x3 tensor_k = body_k->getInvInertiaTensorWorld();

            // the block-diagonal matrix M = diag(Mj,Mk) // see crash_course.pdf 2.2 end, page 14
            cpMatrix M_inv(12,12); // inverse mass matrix
            M_inv.initializeBlock(mI_j, 0, 0);
            M_inv.initializeBlock(tensor_j, 3, 3);
            M_inv.initializeBlock(mI_k, 6, 6);
            M_inv.initializeBlock(tensor_k, 9, 9);

            btScalar S1 = (G1*M_inv*G1.transpose())(0,0); // The Product returns a 1x1 Matrix, which is why we just save S as a scalar
            btScalar S2 = (G2*M_inv*G2.transpose())(0,0);

            //Calculate Impulse
            btVector12 u = {body_j->getLinearVelocity(), body_j->getAngularVelocity(), body_k->getLinearVelocity(), body_k->getAngularVelocity()};
            cpMatrix u_mat(12,1); //12x1
            for(int i = 0; i < 4; i++){
                u_mat.setWithBtVector3(u[i], i*3, 0);
            }


            auto mu = getMU();
            auto lambda = contact.m_appliedImpulse;
            btScalar deltaA1 = 0;
            btScalar deltaA2 = 0;
            if(S1 != 0){
                // deltaA1 = ((G1 * u_mat) * -1 * (1/S1))(0,0);
                //a1 += deltaA1;
                //cout << "before: " << dec << a1 << endl;
                // if(a1+deltaA1 < -0.5*contact.m_appliedImpulse){
                //     deltaA1 = -0.5*contact.m_appliedImpulse - a1;
                //     a1 = -0.5*contact.m_appliedImpulse;
                // }else if(a1+deltaA1 > 0.5*contact.m_appliedImpulse){
                //     deltaA1 = 0.5*contact.m_appliedImpulse - a1;
                //     a1 = 0.5*contact.m_appliedImpulse;
                // }else{
                //     a1 += deltaA1;
                // }
                // //cout << "after: " << dec << a1 << endl;

                auto& a1 = contact.m_appliedImpulseLateral1;
                
                deltaA1 = ((G1 * u_mat) * -1 * (1/S1))(0,0);

                btClamp<btScalar>(deltaA1, -1 * (mu * lambda) - a1, (mu * lambda) - a1);

                //cout << "a1: " << dec << a1 << endl;
                //cout << "delta a1: " << dec << deltaA1 << endl;

                a1 += deltaA1;

                //cout << "mu*lambda: " << dec << mu * lambda << endl;

                //assert(a1 >= -1 * (mu * lambda));
                //assert(a1 <= (mu * lambda));
            }
            if(S2 != 0){        
                auto& a2 = contact.m_appliedImpulseLateral2;

                deltaA2 = ((G2 * u_mat) * -1 * (1/S2))(0,0);

                //cout << "a2: " << dec << a2 << endl;
                //cout << "delta a2: " << dec << deltaA2 << endl;


                btClamp<btScalar>(deltaA2, -1 * (mu * lambda) - a2, (mu * lambda) - a2);

                a2 += deltaA2;
                
                //cout << "mu*lambda: " << dec << mu * lambda << endl;

                //assert(a2 >= -1 * (mu * lambda));
                //assert(a2 <= (mu * lambda));
            }

            cpMatrix M_invG_transposeDeltaLambda1 = M_inv * G1.transpose() * deltaA1; // 12x1
            cpMatrix M_invG_transposeDeltaLambda2 = M_inv * G2.transpose() * deltaA2; // 12x1

            //M_invG_transposeDeltaLambda1.print("M_invG...");
            body_j->setLinearVelocity( u[0] + M_invG_transposeDeltaLambda1.getBtVector3(0,0) + M_invG_transposeDeltaLambda2.getBtVector3(0,0));
            body_j->setAngularVelocity(u[1] + M_invG_transposeDeltaLambda1.getBtVector3(3,0) + M_invG_transposeDeltaLambda2.getBtVector3(3,0));
            body_k->setLinearVelocity( u[2] + M_invG_transposeDeltaLambda1.getBtVector3(6,0) + M_invG_transposeDeltaLambda2.getBtVector3(6,0));
            body_k->setAngularVelocity(u[3] + M_invG_transposeDeltaLambda1.getBtVector3(9,0) + M_invG_transposeDeltaLambda2.getBtVector3(9,0));
        }
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