#ifndef CUSTOM_DYNAMICS_WORLD_H
#define CUSTOM_DYNAMICS_WORLD_H

#include "cpMatrix.h"
#include <BulletSoftBody/btSoftRigidDynamicsWorld.h>
#include <BulletDynamics/ConstraintSolver/btPoint2PointConstraint.h>
#include <BulletDynamics/ConstraintSolver/btHingeConstraint.h>

#include <vector>
#include <iostream>
#include <string>
#include <cassert>

using namespace std;

#define I btMatrix3x3::getIdentity()

#define epsilon 0.0000001f

class CustomDynamicsWorld : public btSoftRigidDynamicsWorld {
private:
    int m_constraint_iters;
    int m_hinge_iters;
    float m_contact_gamma; // attenuator for drift correction for normal correction per iteration
    float m_ball_gamma; // attenuator for drift correction for ball joint correction per iteration
    float m_hinge_gamma; // attenuator for drift correction for hinge joint correction (ball joint part + axis constraint) per iteration
    bool m_apply_friction_constraints; // toggle to correct for friction
    bool m_apply_contact_constraints; // toggle to correct for collisions
    bool m_apply_ball_joints_constraints; // toggle to correct for ball joints
    bool m_apply_hinge_joints_constraints;
    bool m_warm_starting;
    float m_warm_starting_factor;
    bool m_hinge_with_2x2;

public:
    CustomDynamicsWorld(btDispatcher* dispatcher,
                        btBroadphaseInterface* pairCache,
                        btConstraintSolver* constraintSolver,
                        btCollisionConfiguration* collisionConfiguration,
                        btSoftBodySolver* softBodySolver = 0)
                            : btSoftRigidDynamicsWorld(dispatcher, pairCache, constraintSolver, collisionConfiguration, softBodySolver),
                              m_constraint_iters(30), m_hinge_iters(10), m_contact_gamma(0.1f), m_ball_gamma(0.05f), m_hinge_gamma(0.1f), 
                              m_apply_friction_constraints(true), m_apply_ball_joints_constraints(true), m_apply_contact_constraints(true),
                              m_apply_hinge_joints_constraints(true), m_warm_starting(true), m_warm_starting_factor(1.f), m_hinge_with_2x2(false) {}

    void setConstraintIterations(int iterations) { m_constraint_iters = iterations; }
    int getConstraintIterations() const { return m_constraint_iters; }

    void setHingeIterations(int iterations) { m_hinge_iters = iterations; }
    int getHingeIterations() const { return m_hinge_iters; }

    void setContactGamma(float gamma) { m_contact_gamma = gamma; }
    float getContactGamma() const { return m_contact_gamma; }

    void setBallGamma(float gamma) { m_ball_gamma = gamma; }
    float getBallGamma() const { return m_ball_gamma; }

    void setHingeGamma(float gamma) { m_hinge_gamma = gamma; }
    float getHingeGamma() const { return m_hinge_gamma; }

    void setApplyFrictionCorrections(bool friction) { m_apply_friction_constraints = friction; }
    bool getApplyFrictionCorrections() const { return m_apply_friction_constraints; }

    void setApplyContactCorrections(bool contact) { m_apply_contact_constraints = contact; }
    bool getApplyContactCorrections() const { return m_apply_contact_constraints; }

    void setApplyBallJointsCorrections(bool ball_joints) { m_apply_ball_joints_constraints = ball_joints; }
    bool getApplyBallJointsCorrections() const { return m_apply_ball_joints_constraints; }

    void setApplyHingeJointsCorrections(bool hinge_joints) { m_apply_hinge_joints_constraints = hinge_joints; }
    bool getApplyHingeJointsCorrections() const { return m_apply_hinge_joints_constraints; }

    void setWarmStarting(bool warm_starting) { m_warm_starting = warm_starting; }
    bool getWarmStarting() const { return m_warm_starting; }

    void setWarmStartingFactor(float factor) { m_warm_starting_factor = factor; }
    float getWarmStartingFactor() const { return m_warm_starting_factor; }

    void setHingeWith2x2(bool hw2) { m_hinge_with_2x2 = hw2; }
    bool getHingeWith2x2() const { return m_hinge_with_2x2; }

protected:
    void internalSingleStepSimulation(btScalar timeStep) override;

	void sequentialImpulses(btScalar timeStep);

	void point2PointConstraintCorrection(std::vector<btPoint2PointConstraint *> &constraints, btScalar timeStep);

    void hingeJointConstraintCorrection(std::vector<btHingeConstraint *> &constraints, btScalar timeStep, vector<btScalar> &accumulated_impulses);

    void hingeBallJointConstraint(btHingeConstraint* c, btScalar timeStep);
    void hingeIndividualAxisConstraint(btHingeConstraint* c, btScalar timeStep);
    void hingeCombinedAxisConstraint(btHingeConstraint* c, btScalar timeStep);
    void hingeMotorConstraint(btHingeConstraint* c, int c_ind, btScalar timeStep, vector<btScalar> &accumulated_impulses);

    void manifoldCorrection(vector<btPersistentManifold *> &manifolds, btScalar timeStep, int iteration,  vector<btScalar> &target_velocities);

	//
	void integrateConstrainedBodiesWithCustomPhysics(btScalar timeStep);

    // Helper functions:
    btScalar multiplyVector3withMatrix3x3FromBothSides(btVector3& vec, btMatrix3x3& mat);
    cpMatrix computeConstraintVelocityMap(btMatrix3x3 R_j, btMatrix3x3 R_k,btVector3 r_j, btVector3 r_k);
    cpMatrix computeInverseMassMatrix(btRigidBody body_j, btRigidBody body_k);
    cpMatrix computeConstraintVelocityMatrix(btRigidBody body_j, btRigidBody body_k);
    void applyImpulse(btRigidBody& body_j, btRigidBody& body_k, cpMatrix impulse);
    bool isZero(btScalar z){ return (z > -epsilon && z < epsilon);}
    bool isZero(btVector3 v){ return (v.length2() > -epsilon && v.length2() < epsilon);}

    //Prints:
    void printVector(const btVector3& vector, string name);
    void printMatrix(const btMatrix3x3& matrix, string name);
    void printQuat(const btQuaternion quat, string name);
};

#endif // CUSTOM_DYNAMICS_WORLD_H
