#include "custom_dynamics_world.h"

typedef btVector3 btVector12[4];


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

/* sequential impulses method:
    (1. Update velocities of rigid bodies by applying external forces (i.e. set their velocity to u = u^p)) happens earlier
    2. Until convergence or maximum number of iterations:
        2.1 Foreach Constraint i: Compute Si for the system of the two rigid bodies constrained only by the current constraint in isolation
            2.1.1 Compute constraint velocity map G
            2.1.2 getting mass matrix M (or rather M^{-1})
            2.1.3 compute S
        2.2 Compute an impulse represented by ∆λ˜_i by ∆λ˜_i = S^{−1]_i(−Giu), where u = (u_j , u_k).
        2.3 Apply the impulse by updating the velocities of rigid bodies j and k, i.e. u <- u^p + M^{-1}G^T_i ∆λ˜
            (2.3.1 compute u^p = u^n + ∆t  M−1 f^{ext})
            2.3.2 compute M^{-1}G^T_i ∆λ
            2.3.3 compute new velocity u
    (3. Update positions by steps 3.1 and 3.2 (velocities are already up-to-date).
        3.1 Compute new positions z^{n+1} = z^n + ∆t H u^{n+1}  NOTE: (in practice: use quaternion update like before)
        3.2 Normalize quaternions to obtain final rotational state (avoids drift from unit property)
        3.3 Update Transform) */
void CustomDynamicsWorld::sequentialImpulses(std::vector<btPoint2PointConstraint *>& constraints, btScalar timeStep){
    
    // 2
    for (int i = 0; i < getConstraintIterations(); i++){  
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

    // filter out the btPoint2PointConstraint instances:
    const auto num_constraints = this->getNumConstraints();
    std::vector<btPoint2PointConstraint *> point_constraints;

    for (int i = 0; i < num_constraints; ++i){
        auto c = this->getConstraint(i);
        if(c->getConstraintType() == POINT2POINT_CONSTRAINT_TYPE){
            point_constraints.push_back(dynamic_cast<btPoint2PointConstraint *>(c));
        }
    }
 
    sequentialImpulses(point_constraints, timeStep);

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