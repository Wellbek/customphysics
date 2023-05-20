#include "custom_dynamics_world.h"

#include <BulletDynamics/ConstraintSolver/btPoint2PointConstraint.h>

#include <vector>
#include <iostream>

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

void CustomDynamicsWorld::integrateConstrainedBodiesWithCustomPhysics(btScalar timeStep) {
    // We typically perform collision detection at the beginning of the step
    performDiscreteCollisionDetection();

    // Collect all instances of btRigidBody (as pointers) in the current simulation. Note: there may also be
    // collision objects which are not rigid bodies.
    auto bodies = collectBodies(getCollisionObjectArray());

    for (auto body : bodies)
    {
        auto x = body->getCenterOfMassPosition();
        auto q = body->getOrientation();
        body->setLinearVelocity(body->getLinearVelocity() + timeStep * body->getInvMass() * body->getTotalForce());
        // body->setAngularVelocity(body->getAngularVelocity() + body->getInvInertiaTensorWorld());
        x = x + timeStep * body->getLinearVelocity();
        q = q + 1/2 * timeStep * body->getAngularVelocity() * q;
        body->applyGravity();
        body->setCenterOfMassTransform(btTransform(q, x));
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

    // sequential impulses method:
    // 1. Update velocities of rigid bodies by applying external forces. (see above)
    for (auto c : point_constraints){
        btRigidBody* body_j = &c->getRigidBodyA();
        btRigidBody* body_k = &c->getRigidBodyB();
        btMatrix3x3 R_j, R_k; // rotation matrices
        const btVector3* r_j = &c->getPivotInA(); // attachement point for body j
        const btVector3* r_k = &c->getPivotInA(); // attachement point for body k
        const btMatrix3x3 I = btMatrix3x3::getIdentity();
        btVector3 u_j = body_j->getLinearVelocity();
        btVector3 u_k = body_k->getLinearVelocity();
        // 2. Until convergence or maximum number of iterations:
        for (int i = 0; i < getConstraintIterations(); i++){
            // Compute Si for the system of the two rigid bodies constrained only by the current constraint in isolation
            btMatrix3x3 R_j = body_j->getWorldTransform().getBasis();
            btMatrix3x3 R_k = body_k->getWorldTransform().getBasis();

            //auto G = (I, - (R_j * *r_j), I * -1, (R_k * *r_k)); 
            auto G = I; // TODO: replace I by actual G computation

            // NOTE: CHAT GPT und so...
            auto inv_mass_j = body_j->getInvMass();
            auto inv_mass_k = body_k->getInvMass();
            auto inv_total_mass = inv_mass_j + inv_mass_k;
            auto inv_inertia_j = body_j->getInvInertiaTensorWorld();
            auto inv_inertia_k = body_k->getInvInertiaTensorWorld();
            auto inv_total_inertia = inv_inertia_j + inv_inertia_k;
            auto M = inv_total_inertia.inverse() + I * inv_total_mass; // mass matrix

            auto S = G * M.inverse() * G.transpose();
            // Compute an impulse represented by ∆λ˜_i by ∆λ˜_i = S^{−1]_i(−Giu), where u = (u_j , u_k).
            btVector3 impulse_j = S.inverse() * (G * -1 * u_j);
            btVector3 impulse_k = S.inverse() * (G * -1 * u_k);
            // Apply the impulse by updating the velocities of rigid bodies j and k, i.e. u <- u + M^{-1}G^T_i ∆λ˜
            u_j += M.inverse() * G.transpose() * impulse_j;
            u_k += M.inverse() * G.transpose() * impulse_k;
            body_j->setLinearVelocity(u_j);
            body_k->setLinearVelocity(u_k);
        }
        // 3. Update positions by steps 3.1 and 3.2 (velocities are already up-to-date).
        // 3.1 Compute new positions z^{n+1} = z^n + ∆t H u^{n+1}  NOTE: (in practice: use quaternion update like before)
        btVector3 z_j = body_j->getCenterOfMassPosition();
        btVector3 z_k = body_k->getCenterOfMassPosition();
        auto H = I;// TODO: replace I by actual H computation
        z_j += H * u_j * timeStep;
        z_k += H * u_k * timeStep;

        // 3.2 Normalize quaternions to obtain final rotational state (avoids drift from unit property)
        btQuaternion q_j = body_j->getOrientation();
        btQuaternion q_k = body_k->getOrientation();
        q_j.safeNormalize();
        q_k.safeNormalize();

        // finally, update positions
        body_j->setCenterOfMassTransform(btTransform(q_j, z_j));
        body_k->setCenterOfMassTransform(btTransform(q_k, z_k));
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
