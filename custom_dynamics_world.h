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
    float m_apply_friction_constraints; // toggle to correct for friction
    float m_apply_contact_constraints; // toggle to correct for collisions
    float m_apply_ball_joints_constraints; // toggle to correct for ball joints

public:
    CustomDynamicsWorld(btDispatcher* dispatcher,
                        btBroadphaseInterface* pairCache,
                        btConstraintSolver* constraintSolver,
                        btCollisionConfiguration* collisionConfiguration,
                        btSoftBodySolver* softBodySolver = 0)
                            : btSoftRigidDynamicsWorld(dispatcher, pairCache, constraintSolver, collisionConfiguration, softBodySolver),
                              m_constraint_iters(10), m_gamma(0.1f), m_mu(0.3f), 
                              m_apply_friction_constraints(true), m_apply_ball_joints_constraints(true), m_apply_contact_constraints(true) {}

    void setConstraintIterations(int iterations) { m_constraint_iters = iterations; }
    int getConstraintIterations() const { return m_constraint_iters; }

    void setGamma(float gamma) { m_gamma = gamma; }
    float getGamma() const { return m_gamma; }
    
    void setMU(float mu) { m_mu = mu; }
    float getMU() const { return m_mu; }

    void setApplyFrictionCorrections(bool friction) { m_apply_friction_constraints = friction; }
    bool getApplyFrictionCorrections() const { return m_apply_friction_constraints; }

    void setApplyContactCorrections(bool contact) { m_apply_contact_constraints = contact; }
    bool getApplyContactCorrections() const { return m_apply_contact_constraints; }

    void setApplyBallJointsCorrections(bool ball_joints) { m_apply_ball_joints_constraints = ball_joints; }
    bool getApplyBallJointsCorrections() const { return m_apply_ball_joints_constraints; }

protected:
    void internalSingleStepSimulation(btScalar timeStep) override;

	void sequentialImpulses(btScalar timeStep);

	void point2PointConstraintCorrection(std::vector<btPoint2PointConstraint *> &constraints, btScalar timeStep);

    void manifoldCorrection(vector<btPersistentManifold *> &manifolds, btScalar timeStep, int i);

	//
	void integrateConstrainedBodiesWithCustomPhysics(btScalar timeStep);

    // Helper functions:
    btScalar multiplyVector3withMatrix3x3FromBothSides(btVector3& vec, btMatrix3x3& mat);
    cpMatrix computeConstraintVelocityMap(btMatrix3x3 R_j, btMatrix3x3 R_k,btVector3 r_j, btVector3 r_k);
    cpMatrix computeInverseMassMatrix(btRigidBody body_j, btRigidBody body_k);
    cpMatrix computeConstraintVelocityMatrix(btRigidBody body_j, btRigidBody body_k);
    void applyImpulse(btRigidBody& body_j, btRigidBody& body_k, cpMatrix impulse);

    //Prints:
    void printVector(const btVector3& vector, string name);
    void printMatrix(const btMatrix3x3& matrix, string name);
    void printQuat(const btQuaternion quat, string name);
};

# define I btMatrix3x3::getIdentity()

#endif // CUSTOM_DYNAMICS_WORLD_H
