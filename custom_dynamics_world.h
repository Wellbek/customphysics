#ifndef CUSTOM_DYNAMICS_WORLD_H
#define CUSTOM_DYNAMICS_WORLD_H

#include <BulletSoftBody/btSoftRigidDynamicsWorld.h>
#include <BulletDynamics/ConstraintSolver/btPoint2PointConstraint.h>

#include <vector>
#include <iostream>
#include <string>

class CustomDynamicsWorld : public btSoftRigidDynamicsWorld {
private:
    int m_constraint_iters;
    float m_gamma; // attenuator for drift correction per iteration

public:
    CustomDynamicsWorld(btDispatcher* dispatcher,
                        btBroadphaseInterface* pairCache,
                        btConstraintSolver* constraintSolver,
                        btCollisionConfiguration* collisionConfiguration,
                        btSoftBodySolver* softBodySolver = 0)
                            : btSoftRigidDynamicsWorld(dispatcher, pairCache, constraintSolver, collisionConfiguration, softBodySolver),
                              m_constraint_iters(10), m_gamma(1.0f) {}

    void setConstraintIterations(int iterations) { m_constraint_iters = iterations; }
    int getConstraintIterations() const { return m_constraint_iters; }

    void setGamma(float gamma) { m_gamma = gamma; }
    float getGamma() const { return m_gamma; }

protected:
    void internalSingleStepSimulation(btScalar timeStep) override;

	void sequentialImpulses(std::vector<btPoint2PointConstraint *> &constraints, btScalar timeStep);

	//
	void integrateConstrainedBodiesWithCustomPhysics(btScalar timeStep);
};


#endif // CUSTOM_DYNAMICS_WORLD_H
