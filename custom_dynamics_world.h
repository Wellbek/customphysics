#ifndef CUSTOM_DYNAMICS_WORLD_H
#define CUSTOM_DYNAMICS_WORLD_H

#include "cpMatrix.h"
#include <BulletSoftBody/btSoftRigidDynamicsWorld.h>
#include <BulletDynamics/ConstraintSolver/btPoint2PointConstraint.h>

#include <vector>
#include <iostream>
#include <string>
#include <cassert>

using namespace std;

class CustomDynamicsWorld : public btSoftRigidDynamicsWorld {
private:
    int m_constraint_iters;
    float m_gamma; // attenuator for drift correction per iteration
    float m_mu; // attenuator for friction correction per iteration

public:
    CustomDynamicsWorld(btDispatcher* dispatcher,
                        btBroadphaseInterface* pairCache,
                        btConstraintSolver* constraintSolver,
                        btCollisionConfiguration* collisionConfiguration,
                        btSoftBodySolver* softBodySolver = 0)
                            : btSoftRigidDynamicsWorld(dispatcher, pairCache, constraintSolver, collisionConfiguration, softBodySolver),
                              m_constraint_iters(10), m_gamma(0.1f), m_mu(0.3f) {}

    void setConstraintIterations(int iterations) { m_constraint_iters = iterations; }
    int getConstraintIterations() const { return m_constraint_iters; }

    void setGamma(float gamma) { m_gamma = gamma; }
    float getGamma() const { return m_gamma; }
    
    void setMU(float mu) { m_mu = mu; }
    float getMU() const { return m_mu; }

protected:
    void internalSingleStepSimulation(btScalar timeStep) override;

	void sequentialImpulses(btScalar timeStep);

	void point2PointConstraintCorrection(std::vector<btPoint2PointConstraint *> &constraints, btScalar timeStep);

	void contactCorrection(std::vector<btPersistentManifold *> &constraints, btScalar timeStep);

	void frictionCorrection(std::vector<btPersistentManifold *> &constraints, btScalar timeStep);

	//
	void integrateConstrainedBodiesWithCustomPhysics(btScalar timeStep);

    // Helper functions:
    cpMatrix computeConstraintVelocityMap(btMatrix3x3 R_j, btMatrix3x3 R_k,btVector3 r_j, btVector3 r_k);
    cpMatrix computeInverseMassMatrix(btRigidBody body_j, btRigidBody body_k);
    cpMatrix computeConstraintVelocityMatrix(btRigidBody body_j, btRigidBody body_k);
    void applyImpulse(btRigidBody& body_j, btRigidBody& body_k, cpMatrix impulse);
};

# define I btMatrix3x3::getIdentity()

#endif // CUSTOM_DYNAMICS_WORLD_H
