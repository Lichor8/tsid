//
// Copyright (c) 2017 CNRS
//
// This file is part of tsid
// tsid is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
// tsid is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Lesser Public License for more details. You should have
// received a copy of the GNU Lesser General Public License along with
// tsid If not, see
// <http://www.gnu.org/licenses/>.
//

#include "tsid/math/utils.hpp"
#include "tsid/tasks/task-se3-equality.hpp"
#include "tsid/robots/robot-wrapper.hpp"

namespace tsid
{
  namespace tasks
  {
    using namespace math;
    using namespace trajectories;
    using namespace se3;

    TaskSE3Equality::TaskSE3Equality(const std::string & name,
                                     RobotWrapper & robot,
                                     const std::string & frameName):
      TaskMotion(name, robot),
      m_frame_name(frameName),
      m_constraint(name, 6, robot.nv()),
      m_ref(12,6) // (pos size=12, vel size=6)
    {
      assert(m_robot.model().existFrame(frameName));
      m_frame_id = m_robot.model().getFrameId(frameName);

      m_v_ref.setZero();
      m_a_ref.setZero();
      m_M_ref.setIdentity();
      m_wMl.setIdentity();
      m_p_error_vec.setZero(6);
      m_v_error_vec.setZero(6);
      m_p.resize(12);
      m_v.resize(6);
      m_p_ref.resize(12);
      m_v_ref_vec.resize(6);
      m_Kp.setZero(6);
      m_Kd.setZero(6);
      m_a_des.setZero(6);
      m_J.setZero(6, robot.nv());
    }

    int TaskSE3Equality::dim() const
    {
      //return self._mask.sum ()
      return 6; // did I edit this?
    }

    const Vector & TaskSE3Equality::Kp() const { return m_Kp; }

    const Vector & TaskSE3Equality::Kd() const { return m_Kd; }

    void TaskSE3Equality::Kp(ConstRefVector Kp)
    {
      assert(Kp.size()==6);
      m_Kp = Kp;
    }

    void TaskSE3Equality::Kd(ConstRefVector Kd)
    {
      assert(Kd.size()==6);
      m_Kd = Kd;
    }

    void TaskSE3Equality::setReference(TrajectorySample & ref)
    {
      m_ref = ref;
      vectorToSE3(ref.pos, m_M_ref);
      m_v_ref = Motion(ref.vel);
      m_a_ref = Motion(ref.acc);
    }

    const TrajectorySample & TaskSE3Equality::getReference() const
    {
      return m_ref;
    }

    const Vector & TaskSE3Equality::position_error() const
    {
      return m_p_error_vec;
    }

    const Vector & TaskSE3Equality::velocity_error() const
    {
      return m_v_error_vec;
    }

    const Vector & TaskSE3Equality::position() const
    {
      return m_p;
    }

    const Vector & TaskSE3Equality::velocity() const
    {
      return m_v;
    }

    const Vector & TaskSE3Equality::position_ref() const
    {
      return m_p_ref;
    }

    const Vector & TaskSE3Equality::velocity_ref() const
    {
      return m_v_ref_vec;
    }

    const Vector & TaskSE3Equality::getDesiredAcceleration() const
    {
      return m_a_des;
    }

    Vector TaskSE3Equality::getAcceleration(ConstRefVector dv) const
    {
      return m_constraint.matrix()*dv + m_drift.toVector();
    }

    Index TaskSE3Equality::frame_id() const
    {
      return m_frame_id;
    }

    const ConstraintBase & TaskSE3Equality::getConstraint() const
    {
      return m_constraint;
    }

    const ConstraintBase & TaskSE3Equality::compute(const double ,
                                                    ConstRefVector ,
                                                    ConstRefVector ,
                                                    const Data & data)
    {
      SE3 oMi;            // frame position and orientation w.r.t. world frame
      Motion v_frame;     // frame "twist" w.r.t. to the parent joint
      m_robot.framePosition(data, m_frame_id, oMi);
      m_robot.frameVelocity(data, m_frame_id, v_frame);
      m_robot.frameClassicAcceleration(data, m_frame_id, m_drift);

      // Transformation from local to world
      m_wMl.rotation(oMi.rotation());
      std::cout<<"m_wMl:"<<std::endl;
      std::cout<<m_wMl<<std::endl;

      // m_p_error is a Motion data class, which is essentially holds a 6x1 vector
      // in this case errorInSE3 returns a twist error 6x1 vector during time one
      // which is essentially a conversion from twist to position and orientation?
      errorInSE3(oMi, m_M_ref, m_p_error);          // pos err in local frame
      m_v_error = v_frame - m_wMl.actInv(m_v_ref);  // vel err in local frame

      m_p_error_vec = m_p_error.toVector();   // local frame 6x1
      m_v_error_vec = m_v_error.toVector();   // local frame 6x1
      SE3ToVector(m_M_ref, m_p_ref);          // world frame 12x1
      m_v_ref_vec = m_v_ref.toVector();       // world frame 6x1
      SE3ToVector(oMi, m_p);                  // world frame 12x1
      m_v = v_frame.toVector();               // local frame 6x1

      // debug
      PRINT_VECTOR(m_p);
      PRINT_VECTOR(m_p_ref);
      PRINT_VECTOR(m_p_error_vec);
      std::cout<<v_frame<<std::endl;

#ifndef NDEBUG
//      PRINT_VECTOR(v_frame.toVector());
//      PRINT_VECTOR(m_v_ref.toVector());
#endif

      // desired acc in local frame
      m_a_des = - m_Kp.cwiseProduct(m_p_error.toVector())
                - m_Kd.cwiseProduct(m_v_error.toVector())
                + m_wMl.actInv(m_a_ref).toVector();

      // @todo Since Jacobian computation is cheaper in world frame
      // we could do all computations in world frame
      m_robot.frameJacobianLocal(data, m_frame_id, m_J);

//______________________________________________________________________________
      // task point (C local frame)
      // todo: check order of error
      // todo: make nv variable for DoF size (i.e. handle variable DoF)
      // todo: check drift, is 0?

      // input: m_ref.pos = A^o_E
      // input: m_ref.vel = C[A]^v_A,E
      // input: m_ref.acc = C[A]^vd_A,E

//      TrajectorySample ref;
//      ref.pos = m_ref.pos.head(3);  // get first 3 entries from 12x1 vector world frame
//      ref.vel = m_ref.vel.head(3);  // get first 3 entries from 6x1 vector world frame

//      Motion v_ref = Motion(ref.vel);

//      Vector3 p_error_vec;
//      p_error_vec = oMi.translation() - ref.pos;    // 3x1 - 3x1  pos err in world frame
//      p_error_vec = m_wMl.actInv(p_error_vec);      // 3x1 * 3x3 pos err in local frame
//      Vector v = m_v.head(3);                       // get first 3 entries from 6x1 vector local frame
//      Vector v_error_vec = v - m_wMl.actInv(v_ref).toVector().head(3);  // 3x1 - 3x1 vel err in local frame
//      Vector a_ref = m_wMl.actInv(m_a_ref).toVector().head(3);

//      // desired acc in local frame
//      Vector a_des = - m_Kp.cwiseProduct(p_error_vec)
//                     - m_Kd.cwiseProduct(v_error_vec)
//                     + a_ref;                              // 3x1 desired acc in local frame

//      // use only the 3D part
//      Matrix3x J = m_J.block(0,0,3,3);             // J (3x3)
//      Vector6 drift;
//      drift.head(3) = m_drift.toVector().head(3);  // Jd*qd (3x1)

//      // convert to 6D for se3-task (required)
//      const int nv = 6;
//      m_J.setZero(6, nv);
//      m_a_des.setZero(6);

//      m_J.block(0,0,3,3) = J;
//      m_a_des.head(3) = a_des;

//      // debug
//      std::cout<<"m_J:"<<std::endl<<m_J<<std::endl;
//      std::cout<<"m_a_des:"<<std::endl<<m_a_des<<std::endl;
//      std::cout<<"drift:"<<std::endl<<drift<<std::endl;

//      // in local frame:
//      // ||A*qdd - b||^2
//      m_constraint.setMatrix(m_J); // A (6x6)
//      m_constraint.setVector(m_a_des - drift); // b (6x1)
//______________________________________________________________________________

//______________________________________________________________________________
      // task point (C[A] mixed frame)
      // todo: check order of error
      // todo: make nv variable for DoF size (i.e. handle variable DoF) (use Jacobian size? or just set zero?)
      // todo: check drift, is 0?
      // todo: drift should also use a 3x6 J, however I do not control this part of the software, so
      // minor influence of rotation part? or does this part end up being zero due to other zeor parts?

      // input: m_ref.pos = A^o_E
      // input: m_ref.vel = C[A]^v_A,E
      // input: m_ref.acc = C[A]^vd_A,E

      std::cout<<"--------Task Point-------"<<std::endl;

      m_J = m_wMl.toActionMatrix()*m_J; // Jacobian: C[A]^J_A,C       = C[A]^X_C * C^J_A,C
      m_drift = m_wMl.act(m_drift);     // drift:    C[A]^Jd_A,C * qd = C[A]^X_C * C^Jd_A,C * qd

      TrajectorySample ref;
      ref.pos = m_ref.pos.head(3);  // get first 3 entries from 12x1 vector world frame
      ref.vel = m_ref.vel.head(3);  // get first 3 entries from 6x1 vector world frame
      ref.acc = m_a_ref.toVector().head(3);

      Vector3 p_error_vec;
      p_error_vec = oMi.translation() - ref.pos;    // 3x1 - 3x1  pos err in world frame
      m_v = m_wMl.act(v_frame).toVector();          // m_v required for velocity() method
      Vector3 v_error_vec = m_v.head(3) - ref.vel;  // 3x1 - 3x1 vel err in local[world] frame
      Vector3 a_ref = ref.acc;

      std::cout<<"m_v:"<<std::endl<<m_v<<std::endl;
      std::cout<<"ref.vel:"<<std::endl<<ref.vel<<std::endl;

      std::cout<<"m_a_ref:"<<std::endl<<m_ref.acc<<std::endl;



////       separated task-point hack
////      Motion a_frame;
////      m_robot.frameAcceleration(data, m_frame_id, a_frame);

//      v_error_vec = v_error_vec - m_wMl.act(v_frame).toVector().head(3);
////      a_ref = a_ref + m_wMl.act(a_frame).toVector().head(3);
//      std::cout<<"data.ddq:"<<std::endl<<data.ddq<<std::endl;
//      a_ref = a_ref + (m_J*data.ddq).head(3) + m_drift.toVector().head(3);

      // desired acc in world frame
      Vector3 a_des = - m_Kp.cwiseProduct(p_error_vec)
                     - m_Kd.cwiseProduct(v_error_vec)
                     + a_ref;                              // 3x1 desired acc in local frame

      // use only the 3D part
      Matrix3x J = m_J.block(0,0,3,6);             // J (3x3)
      Vector6 drift;
      drift.head(3) = m_drift.toVector().head(3);  // Jd*qd (3x1)

      // convert to 6D for se3-task (required)
      const int nv = 6;
      m_J.setZero(6, nv);
      m_a_des.setZero(6);

      m_J.block(0,0,3,6) = J;
      m_a_des.head(3) = a_des;

      // debug
      std::cout<<"m_J:"<<std::endl<<m_J<<std::endl;
      std::cout<<"m_a_des:"<<std::endl<<m_a_des<<std::endl;
      std::cout<<"drift:"<<std::endl<<drift<<std::endl;

      PRINT_VECTOR(ref.pos);
      PRINT_VECTOR(oMi.translation());
      PRINT_VECTOR(p_error_vec);

      std::cout<<"-----------------------"<<std::endl;

      // in local frame:
      // ||A*qdd - b||^2
      m_constraint.setMatrix(m_J); // A (6x6)
      m_constraint.setVector(m_a_des - drift); // b (6x1)
//______________________________________________________________________________

//      // in local frame:
//      // ||A*qdd - b||^2
//      m_constraint.setMatrix(m_J); // A (6x6)
//      m_constraint.setVector(m_a_des - m_drift.toVector()); // b (6x1)
      return m_constraint;
    }
    
  }
}
