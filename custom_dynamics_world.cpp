#include "custom_dynamics_world.h"

#include <BulletDynamics/ConstraintSolver/btPoint2PointConstraint.h>

#include <vector>
#include <iostream>
#include <string>

typedef btMatrix3x3 btMatrix12x12[4][4];
typedef btMatrix3x3 btMatrix3x12[4];
typedef btMatrix3x3 btMatrix12x3[4];
typedef btVector3 btVector12[4];

// for a diagonal matrix consisting of blocks with each block being invertible, the inverse of the matrix is given by the inverse of each block matrix
// https://en.wikipedia.org/wiki/Block_matrix
void inverseDiagBlock3x3(btMatrix12x12& mat){
    for (int i = 0; i < 4; i++) {
        btMatrix3x3& block = mat[i][i];
        block = block.inverse();
    }
}

// transpose of each block + rearrangement of blocks (btMatrix3x12 -> btMatrix12x3)
void transposeConstraintMap(const btMatrix3x12& mat, btMatrix12x3& res){
    for (int i = 0; i < 4; i++) {
        res[i] = mat[i].transpose();
    }
}

void multiply(const btMatrix3x12& mat1, const btMatrix12x3& mat2, btMatrix3x3& res){
    int i,j,x;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            res[i][j] = 0;
            for (x = 0; x < 12; x++) {
                res[i][j] += mat1[x%4][i][x%3] * mat2[x%4][x%3][j];
            }
        }
    }
}

void multiply(const btMatrix12x3& mat1, const btMatrix3x3& mat2, btMatrix12x3& res){
    int i,j,x;
    for (i = 0; i < 12; i++) {
        for (j = 0; j < 3; j++) {
            res[i%4][i%3][j] = 0;
            for (x = 0; x < 3; x++) {
                res[i%4][i%3][j] += mat1[i%4][i%3][x] * mat2[x][j];
            }
        }
    }
}

void multiply(const btMatrix3x3& mat1, const btMatrix3x12& mat2, btMatrix3x12& res){
    int i,j,x;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 12; j++) {
            res[j%4][i][j%3] = 0;
            for (x = 0; x < 3; x++) {
                res[j%4][i][j%3] += mat1[i][x] * mat2[j%4][x][j%3];
            }
        }
    }
}

void multiply(const btMatrix3x12& mat1, const btMatrix12x12& mat2, btMatrix12x3& res){
    int i,j,x;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 12; j++) {
            res[j%4][i][j%3] = 0;
            for (x = 0; x < 12; x++) {
                res[j%4][i][j%3] += mat1[x%4][i][x%3] * mat2[x%4][j%4][x%3][j%3];
            }
        }
    }
}

void multiply(const btMatrix12x12& mat1, const btMatrix12x3& mat2, btMatrix3x12& res){
    int i,j,x;
    for (i = 0; i < 12; i++) {
        for (j = 0; j < 3; j++) {
            res[i%4][i%3][j] = 0;
            for (x = 0; x < 12; x++) {
                res[i%4][i%3][j] += mat1[i%4][x%4][i%3][x%3] * mat2[x%4][x%3][j];
            }
        }
    }
}

void multiply(const btMatrix3x12& mat, btScalar s, btMatrix3x12& res){
    int i,j;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 12; j++) {
            res[j%4][i][j%3] = mat[j%4][i][j%3]*s;
        }
    }
}

void multiply(const btMatrix3x12& mat, const btVector12& vec, btVector3& res){
    int i,x;
    for (i = 0; i < 3; i++) {
        res[i] = 0;
        for (x = 0; x < 12; x++) {
            res[i] += mat[x%4][i][x%3] * vec[x%4][x%3];
        }
    }
}

void multiply(const btMatrix12x3& mat, const btVector3& vec, btVector12& res){
    int i,x;
    for (i = 0; i < 12; i++) {
        res[i%4][i%3] = 0;
        for (x = 0; x < 3; x++) {
            res[i%4][i%3] += mat[i%4][i%3][x] * vec[x];
        }
    }
}

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


void printMatrix(const btMatrix3x12& matrix, char name) {
    std::cout << name << ":"; 
    for(int i = 0; i < 4; i++){
        printMatrix(matrix[i], ' ');
    }
}

void printMatrix(const btMatrix12x12& matrix, char name) {
    std::cout << name << ":"; 
    for (int row = 0; row < 4; row++) {
        for (int col = 0; col < 4; col++) {
            printMatrix(matrix[row][col], ' ');
        }
    }
}
/*
void printMatrix(const btMatrix12x3& matrix, char name) {
    std::cout << name << ":"; 
    for (int row = 0; row < 4; row++) {
        printMatrix(matrix[row], ' ');
    }
}*/

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
        btRigidBody body_j = c->getRigidBodyA();
        btRigidBody body_k = c->getRigidBodyB();
        btMatrix3x3 R_j, R_k; // rotation matrices
        const btVector3* r_j = &c->getPivotInA(); // attachement point for body j
        const btVector3* r_k = &c->getPivotInB(); // attachement point for body k
        const btMatrix3x3 I = btMatrix3x3::getIdentity();
        // 2. Until convergence or maximum number of iterations:
        for (int i = 0; i < getConstraintIterations(); i++){
            // 2.1 Compute Si for the system of the two rigid bodies constrained only by the current constraint in isolation
            //  2.1.1 Compute constraint velocity map G
            //  2.1.2 getting mass matrix M (or rather M^{-1})
            //  2.1.3 compute S
            btMatrix3x3 R_j = (body_j).getWorldTransform().getBasis();
            btMatrix3x3 R_k = (body_k).getWorldTransform().getBasis();
            //printMatrix(R_j);

            // 2.1.1
            btVector3 v1,v2,v3;
            (R_j * *r_j).getSkewSymmetricMatrix(&v1,&v2,&v3);
            //printVector(*r_j);
            //printVector(R_j * *r_j);
            //printVector(v1, 'x'); printVector(v2, 'y'); printVector(v3, 'z');
            btMatrix3x3 K_j = btMatrix3x3(v1,v2,v3);
            //printMatrix(K_j, 'K');
            (R_k * *r_k).getSkewSymmetricMatrix(&v1,&v2,&v3);
            btMatrix3x3 K_k = btMatrix3x3(v1,v2,v3);
            //btMatrix3x12 G = {I,  K_j * -1, I * -1, K_k}; 
            btMatrix3x12 G;
            G[0] = I;
            G[1] = K_j * -1;
            G[2] = I * -1;
            G[3] = K_j;
            //printMatrix(G, 'G');
            btMatrix12x3 G_transpose;
            transposeConstraintMap(G, G_transpose);
            //printMatrix(G_transpose, 'T');
        
            // 2.1.2
            btMatrix3x3 e; // null matrix
            btMatrix3x3 mI_j = I * body_j.getInvMass();
            btMatrix3x3 mI_k = I * body_k.getInvMass();
            btMatrix3x3 tensor_j = body_j.getInvInertiaTensorWorld();
            btMatrix3x3 tensor_k = body_k.getInvInertiaTensorWorld();
            // the block-diagonal matrix M = diag(Mj,Mk) // see crash_course.pdf 2.2 end, page 14
            btMatrix12x12 M_inv = {
                { mI_j, e, e, e},
                { e, tensor_j, e, e},
                { e, e, mI_k, e},
                { e, e, e, tensor_k}
            }; // inverse mass matrix
            //printMatrix(M_inv, 'M');

            // 2.1.3           
            btMatrix3x12 GM_inv;
            multiply(G, M_inv, GM_inv);
            btMatrix3x3 S;
            multiply(GM_inv, G_transpose, S);
            //printMatrix(S, 'S');
            //printMatrix(GM_inv, '3');

            // 2.2 Compute an impulse represented by ∆λ˜_i by ∆λ˜_i = S^{−1]_i(−Giu), where u = (u_j , u_k).
            btVector12 u = {body_j.getLinearVelocity(), body_j.getAngularVelocity(), body_k.getLinearVelocity(), body_k.getAngularVelocity()};
            //printVector(u[0], '1'); printVector(u[1], '2'); printVector(u[2], '3'); printVector(u[3], '4');
            btMatrix3x12 negG;
            multiply(G, -1, negG);
            //printMatrix(negG, '-');
            btVector3 ndC;
            multiply(negG, u, ndC);
            //printVector(ndC, 'C');
            btVector3 impulse;
            if(S.determinant() != 0){
                impulse = S.inverse() * ndC;
            }else{
                impulse = {0, 0, 0};
            }
            //printVector(impulse, 'i');
            // 2.3 Apply the impulse by updating the velocities of rigid bodies j and k, i.e. u <- u + M^{-1}G^T_i ∆λ˜
            btMatrix12x3 M_invG_transpose;
            multiply(M_inv, G_transpose, M_invG_transpose);
            multiply(M_invG_transpose, impulse, u);

            c->getRigidBodyA().setLinearVelocity(body_j.getLinearVelocity() + u[0]);
            c->getRigidBodyA().setAngularVelocity(body_j.getAngularVelocity() + u[1]);
            c->getRigidBodyB().setLinearVelocity(body_k.getLinearVelocity() + u[2]);
            c->getRigidBodyB().setAngularVelocity(body_k.getAngularVelocity() + u[3]);
        }
        // 3. Update positions by steps 3.1 and 3.2 (velocities are already up-to-date).
        // 3.1 Compute new positions z^{n+1} = z^n + ∆t H u^{n+1}  NOTE: (in practice: use quaternion update like before)
        btVector3 x_j = body_j.getCenterOfMassPosition();
        x_j += timeStep * body_j.getLinearVelocity();
        btVector3 x_k = body_k.getCenterOfMassPosition();
        x_k += timeStep * body_k.getLinearVelocity();

        // TODO update rotation
        btVector3 omega_j = body_j.getAngularVelocity();
        btQuaternion q_j = body_j.getOrientation();
        //q_j += timeStep *  btQuaternion(omega_j.x(), omega_j.y(), omega_j.z(), 0);
        btVector3 omega_k = body_k.getAngularVelocity();
        btQuaternion q_k = body_k.getOrientation();
        //q_k += timeStep *  btQuaternion(omega_k.x(), omega_k.y(), omega_k.z(), 0);

        // 3.2 Normalize quaternions to obtain final rotational state (avoids drift from unit property)
        q_j.safeNormalize();
        q_k.safeNormalize();

        // finally, update positions
        body_j.setCenterOfMassTransform(btTransform(q_j, x_j));
        body_k.setCenterOfMassTransform(btTransform(q_k, x_k));
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