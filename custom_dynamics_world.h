#ifndef CUSTOM_DYNAMICS_WORLD_H
#define CUSTOM_DYNAMICS_WORLD_H

#include <BulletSoftBody/btSoftRigidDynamicsWorld.h>

class CustomDynamicsWorld : public btSoftRigidDynamicsWorld {
private:
    int m_constraint_iters;

public:
    CustomDynamicsWorld(btDispatcher* dispatcher,
                        btBroadphaseInterface* pairCache,
                        btConstraintSolver* constraintSolver,
                        btCollisionConfiguration* collisionConfiguration,
                        btSoftBodySolver* softBodySolver = 0)
                            : btSoftRigidDynamicsWorld(dispatcher, pairCache, constraintSolver, collisionConfiguration, softBodySolver),
                              m_constraint_iters(10) {}

    void setConstraintIterations(int iterations) { m_constraint_iters = iterations; }
    int getConstraintIterations() const { return m_constraint_iters; }

protected:
    void internalSingleStepSimulation(btScalar timeStep) override;

    //
    void integrateConstrainedBodiesWithCustomPhysics(btScalar timeStep);
};


#endif // CUSTOM_DYNAMICS_WORLD_H
