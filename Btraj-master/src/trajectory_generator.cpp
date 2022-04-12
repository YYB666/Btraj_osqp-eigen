#include "trajectory_generator.h"
using namespace std;    
using namespace Eigen;



int TrajectoryGenerator::BezierPloyCoeffGeneration(
            const vector<Cube> &corridor,
            const MatrixXd &MQM,
            const MatrixXd &pos,
            const MatrixXd &vel,
            const MatrixXd &acc,
            const double maxVel,
            const double maxAcc,
            const int traj_order,
            const double minimize_order,
            const double margin,
            const bool & isLimitVel,
            const bool & isLimitAcc,
            double & obj,
            MatrixXd & PolyCoeff)  // define the order to which we minimize.   1 -- velocity, 2 -- acceleration, 3 -- jerk, 4 -- snap  
{   
#define ENFORCE_VEL  isLimitVel // whether or not adding extra constraints for ensuring the velocity feasibility
#define ENFORCE_ACC  isLimitAcc // whether or not adding extra constraints for ensuring the acceleration feasibility

    double initScale = corridor.front().t; //获取第一段轨迹时间
    double lstScale  = corridor.back().t;
    int segment_num  = corridor.size();

    int n_poly = traj_order + 1;
    int s1d1CtrlP_num = n_poly;  //一段轨迹在一个方向的控制点个数
    int s1CtrlP_num   = 3 * s1d1CtrlP_num; //一段轨迹的总控制点个数

    int equ_con_s_num = 3 * 3; // p, v, a in x, y, z axis at the start point
    int equ_con_e_num = 3 * 3; // p, v, a in x, y, z axis at the end point
    int equ_con_continuity_num = 3 * 3 * (segment_num - 1);  //中间点 p v a 连续性约束
    int equ_con_num   = equ_con_s_num + equ_con_e_num + equ_con_continuity_num; // p, v, a in x, y, z axis in each segment's joint position
    
    int vel_con_num = 3 *  traj_order * segment_num; //速度 控制点个数
    int acc_con_num = 3 * (traj_order - 1) * segment_num;  //加速度 控制点个数

    if( !ENFORCE_VEL )  
        vel_con_num = 0;

    if( !ENFORCE_ACC )
        acc_con_num = 0;

    int high_order_con_num = vel_con_num + acc_con_num; 
    //int high_order_con_num = 0; //3 * traj_order * segment_num;

    int con_num   = equ_con_num + high_order_con_num;  //总约束个数
    int ctrlP_num = segment_num * s1CtrlP_num;  //所有段轨迹的控制点总数

    
    
    /* con_bdk 存储每个约束的类型及约束值 */
  
    vector< pair<double, double> > con_bdk;

    
    if(ENFORCE_VEL)
    {
        /***  Stack the bounding value for the linear inequality for the velocity constraints  ***/
        for(int i = 0; i < vel_con_num; i++)
        {
            pair<double, double>  cb_ie = make_pair( - maxVel,  + maxVel) ;
            con_bdk.push_back(cb_ie);   
        }
    }

    if(ENFORCE_ACC)
    {
        /***  Stack the bounding value for the linear inequality for the acceleration constraints  ***/
        for(int i = 0; i < acc_con_num; i++)
        {
            pair<double, double>  cb_ie = make_pair( - maxAcc,  maxAcc) ; 
            con_bdk.push_back(cb_ie);   
        }
    }

    /*等式约束的顺序
    [s_p_x, s_p_y, s_p_z, s_v_x, s_v_y, s_v_z, s_a_x, s_a_y, s_a_z, e_p_x, e_p_y, e_p_z, e_v_x, e_v_y, e_v_z, e_a_x, e_a_y, e_a_z, 0, 0, ··· ··· 0, 0}'
     起点 位置 速度 加速度 等式约束
     终点 位置 速度 加速度 等式约束
     中间点 位置 速度 加速度 等式约束
    */

    ROS_WARN("[Bezier Trajectory] equality bound %d", equ_con_num);
    for(int i = 0; i < equ_con_num; i ++ ){ 
        double beq_i;
        if(i < 3)                    beq_i = pos(0, i);      //start point: x y z
        else if (i >= 3  && i < 6  ) beq_i = vel(0, i - 3);  //start point: v_x v_y v_z
        else if (i >= 6  && i < 9  ) beq_i = acc(0, i - 6);  //start point: a_x a_y a_z
        else if (i >= 9  && i < 12 ) beq_i = pos(1, i - 9 ); //end point: x y z
        else if (i >= 12 && i < 15 ) beq_i = vel(1, i - 12); //end point: v_x v_y v_z
        else if (i >= 15 && i < 18 ) beq_i = acc(1, i - 15); //end point: a_x a_y a_z
        else beq_i = 0.0;

        pair<double, double> cb_eq = make_pair( beq_i, beq_i ) ; // # cb_eq means: constriants boundary of equality constrain
        con_bdk.push_back(cb_eq); 
    }
    
  
    /* ## define a container for control points' boundary and boundkey ## */ 
    /* ## dataType in one tuple is : boundary type, lower bound, upper bound ## */
    vector< pair<double, double> > var_bdk; 

    for(int k = 0; k < segment_num; k++)  // segment_num 段轨迹
    {   
        Cube cube_     = corridor[k];  //获取每段轨迹所对应的走廊参数
        double scale_k = cube_.t;    //scale factor 将每段轨迹时间映射到 [ 0,1] 区间

        for(int i = 0; i < 3; i++ )  //三个方向 x y z 
        {   
            for(int j = 0; j < n_poly; j ++ ) // n_poly 个控制点
            {   
                

                double lo_bound, up_bound;
                
                if(k > 0)
                {
                    lo_bound = (cube_.box[i].first  + margin) / scale_k;  // margin 立方体边距 
                    up_bound = (cube_.box[i].second - margin) / scale_k;  //对应论文 公式（5） scale_k 本应在等式左侧，此处在等式两端同时除以 scale_k
                }
                else //第一个cube是基准
                {
                    lo_bound = (cube_.box[i].first)  / scale_k;
                    up_bound = (cube_.box[i].second) / scale_k;
                }

                pair<double, double> vb_x  = make_pair( lo_bound, up_bound ) ; // # vb_x means: varialbles boundary of unknowns x (Polynomial coeff)

                var_bdk.push_back(vb_x);
            }
        } 
    }
    
    
    /************************************** OSQP-EIGEN *************************************/
    int numberOfVariables=ctrlP_num;
    int numberOfConstrains=ctrlP_num+con_num;
    // int numberOfConstrains=con_bdk.size()+var_bdk.size();   

    // allocate QP problem matrices and vectores
    Eigen::SparseMatrix<double> hessian;
    Eigen::VectorXd gradient;
    Eigen::SparseMatrix<double> linearMatrix;
    Eigen::VectorXd lowerBound;
    Eigen::VectorXd upperBound;
    
    
    gradient.resize(numberOfVariables);
    for(int m=0;m<numberOfVariables;m++)
    {
        gradient[m]=0.0;
    }
    

    lowerBound.resize(numberOfConstrains);
    upperBound.resize(numberOfConstrains);
    
    // ROS_WARN("set variables boundary");
    
    
    // Set the bounds on constraints.  此处对应 论文 公式（11）
    //   for i=1, ...,con_num : blc[i] <= constraint i <= buc[i] 
    int idx=0;
    
    for(int i = 0; i < con_num; i++ ) {
        
        lowerBound(idx)=con_bdk[i].first;   // Numerical value of lower bound.
        upperBound(idx)=con_bdk[i].second;   // Numerical value of upper bound.
        idx++;
    }
    
    
    for(int j = 0; j<ctrlP_num ; ++j){ 
         
        
        lowerBound(idx)=var_bdk[j].first;   // Numerical value of lower bound.
        upperBound(idx)=var_bdk[j].second;   // Numerical value of upper bound.
        idx++;        
    } 
    
    
    /***************************  Matrix A *************************************/
    /**************************   AC<=B  ***************************************/

    linearMatrix.resize(numberOfConstrains,numberOfVariables);
   
    
    // ROS_WARN("[Bezier Trajectory] Start stacking the Linear Matrix A, inequality part");
    int row_idx = 0; //Matrix A 的行号
    // The velocity constraints
    if(ENFORCE_VEL)
    {   
        for(int k = 0; k < segment_num ; k ++ )
        {   
            for(int i = 0; i < 3; i++)
            {  // for x, y, z loop
                for(int p = 0; p < traj_order; p++)
                {
                    int nzi = 2;
                    int asub[nzi];
                    double aval[nzi];

                    aval[0] = -1.0 * traj_order;  //此处对应 论文 公式（12）
                    aval[1] =  1.0 * traj_order;

                    asub[0] = k * s1CtrlP_num + i * s1d1CtrlP_num + p;    // 此处对应约束在MOSEK中的序号
                    asub[1] = k * s1CtrlP_num + i * s1d1CtrlP_num + p + 1;    //s1CtrlP_num :一段轨迹的控制点数 = 3* s1d1CtrlP_num：一段轨迹在一个方向的控制点数 
                    
                    
                    for(int j=0;j<nzi;j++)
                    {
                        linearMatrix.insert(row_idx,asub[j])=aval[j];
                    }
                    
        
                    row_idx ++;
                }
            }
        }
    }

    // The acceleration constraints
    if(ENFORCE_ACC)
    {
        for(int k = 0; k < segment_num ; k ++ )
        {
            for(int i = 0; i < 3; i++)
            { 
                for(int p = 0; p < traj_order - 1; p++)
                {    
                    int nzi = 3;
                    int asub[nzi];
                    double aval[nzi];

                    aval[0] =  1.0 * traj_order * (traj_order - 1) / corridor[k].t;  //此处对应 论文 公式（12）
                    aval[1] = -2.0 * traj_order * (traj_order - 1) / corridor[k].t;
                    aval[2] =  1.0 * traj_order * (traj_order - 1) / corridor[k].t;
                    asub[0] = k * s1CtrlP_num + i * s1d1CtrlP_num + p;    
                    asub[1] = k * s1CtrlP_num + i * s1d1CtrlP_num + p + 1;    
                    asub[2] = k * s1CtrlP_num + i * s1d1CtrlP_num + p + 2;    
                    
                    for(int j=0;j<nzi;j++)
                    {
                        linearMatrix.insert(row_idx,asub[j])=aval[j];
                    }
                   
                    row_idx ++;
                }
            }
        }
    }
    /*   Start position  */
    {
        // position :
        for(int i = 0; i < 3; i++)
        {  // loop for x, y, z       
            int nzi = 1;
            int asub[nzi];
            double aval[nzi];
            aval[0] = 1.0 * initScale;
            asub[0] = i * s1d1CtrlP_num;

            for(int j=0;j<nzi;j++)
           {
             linearMatrix.insert(row_idx,asub[j])=aval[j];
           }    
            row_idx ++;
        }
        // velocity :
        for(int i = 0; i < 3; i++)
        {  // loop for x, y, z       
            int nzi = 2;
            int asub[nzi];
            double aval[nzi];
            aval[0] = - 1.0 * traj_order;
            aval[1] =   1.0 * traj_order;
            asub[0] = i * s1d1CtrlP_num;
            asub[1] = i * s1d1CtrlP_num + 1;
            
            for(int j=0;j<nzi;j++)
           {
             linearMatrix.insert(row_idx,asub[j])=aval[j];
           }    
            
            row_idx ++;
        }
        // acceleration : 
        for(int i = 0; i < 3; i++)
        {  // loop for x, y, z       
            int nzi = 3;
            int asub[nzi];
            double aval[nzi];
            aval[0] =   1.0 * traj_order * (traj_order - 1) / initScale;
            aval[1] = - 2.0 * traj_order * (traj_order - 1) / initScale;
            aval[2] =   1.0 * traj_order * (traj_order - 1) / initScale;
            asub[0] = i * s1d1CtrlP_num;
            asub[1] = i * s1d1CtrlP_num + 1;
            asub[2] = i * s1d1CtrlP_num + 2;
            
            for(int j=0;j<nzi;j++)
           {
             linearMatrix.insert(row_idx,asub[j])=aval[j];
           }        
            row_idx ++;
        }
    }      

    /*   End position  */
    // ROS_WARN(" end position");
    {   
        // position :
        for(int i = 0; i < 3; i++)
        {  // loop for x, y, z       
            int nzi = 1;
            int asub[nzi];
            double aval[nzi];
            asub[0] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num;
            aval[0] = 1.0 * lstScale;
            
            for(int j=0;j<nzi;j++)
           {
             linearMatrix.insert(row_idx,asub[j])=aval[j];
           }        
            row_idx ++;
        }
        // velocity :
        for(int i = 0; i < 3; i++)
        { 
            int nzi = 2;
            int asub[nzi];
            double aval[nzi];
            asub[0] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num - 1;
            asub[1] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num;
            aval[0] = - 1.0;
            aval[1] =   1.0;
            
            for(int j=0;j<nzi;j++)
           {
             linearMatrix.insert(row_idx,asub[j])=aval[j];
           }        
            row_idx ++;
        }
        // acceleration : 
        for(int i = 0; i < 3; i++)
        { 
            int nzi = 3;
            int asub[nzi];
            double aval[nzi];
            asub[0] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num - 2;
            asub[1] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num - 1;
            asub[2] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num;
            aval[0] =   1.0 / lstScale;
            aval[1] = - 2.0 / lstScale;
            aval[2] =   1.0 / lstScale;

            for(int j=0;j<nzi;j++)
           {
             linearMatrix.insert(row_idx,asub[j])=aval[j];
           }      
            row_idx ++;
        }
    }

    /*   joint points  */
    /********************************************** 对应论文 公式（10）**********************************************/
    // ROS_WARN(" joint position");
    {
        int sub_shift = 0;
        double val0, val1;
        for(int k = 0; k < (segment_num - 1); k ++ )
        {   
            double scale_k = corridor[k].t;
            double scale_n = corridor[k+1].t;
            // position :
            val0 = scale_k;
            val1 = scale_n;
            for(int i = 0; i < 3; i++)  
            {  // loop for x, y, z
                int nzi = 2;
                int asub[nzi];
                double aval[nzi];

                // This segment's last control point
                aval[0] = 1.0 * val0;
                asub[0] = sub_shift + (i+1) * s1d1CtrlP_num - 1;

                // Next segment's first control point
                aval[1] = -1.0 * val1;
                asub[1] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num;
                
                for(int j=0;j<nzi;j++)
                {
                    linearMatrix.insert(row_idx,asub[j])=aval[j];
                }       
                row_idx ++;
            }
            // velocity
            // val0 = pow(scale_k,1-1)=1
            // val1 = pow(scale_n,1-1)=1
            for(int i = 0; i < 3; i++)
            {  
                int nzi = 4;
                int asub[nzi];
                double aval[nzi];
                
                // This segment's last velocity control point
                aval[0] = -1.0;  
                aval[1] =  1.0;
                asub[0] = sub_shift + (i+1) * s1d1CtrlP_num - 2;    
                asub[1] = sub_shift + (i+1) * s1d1CtrlP_num - 1;   
                // Next segment's first velocity control point
                aval[2] =  1.0;  
                aval[3] = -1.0;

                asub[2] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num;    
                asub[3] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num + 1;

                for(int j=0;j<nzi;j++)
                {
                    linearMatrix.insert(row_idx,asub[j])=aval[j];
                }          
                row_idx ++;
            }
            // acceleration :
            val0 = 1.0 / scale_k;
            val1 = 1.0 / scale_n;
            for(int i = 0; i < 3; i++)
            {  
                int nzi = 6;
                int asub[nzi];
                double aval[nzi];
                
                // This segment's last velocity control point
                aval[0] =  1.0  * val0;
                aval[1] = -2.0  * val0;
                aval[2] =  1.0  * val0;
                asub[0] = sub_shift + (i+1) * s1d1CtrlP_num - 3;    
                asub[1] = sub_shift + (i+1) * s1d1CtrlP_num - 2;   
                asub[2] = sub_shift + (i+1) * s1d1CtrlP_num - 1;   
                // Next segment's first velocity control point
                aval[3] =  -1.0  * val1;
                aval[4] =   2.0  * val1;
                aval[5] =  -1.0  * val1;
                asub[3] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num;    
                asub[4] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num + 1;
                asub[5] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num + 2;

                
                for(int j=0;j<nzi;j++)
                {
                    linearMatrix.insert(row_idx,asub[j])=aval[j];
                }          
                row_idx ++;
            }

            sub_shift += s1CtrlP_num;
        }
    }

    /****************************************** 在矩阵A中增加对 控制变量的约束 此处是MOSEK与OSQP-EIGEN的最大区别之一！！！   **********************************************/
    
    for(int i=0;i<numberOfVariables;i++)
    {
        linearMatrix.insert(row_idx,i)=1;
        
        row_idx++;

    }
    /****************************************   *********************************************************/

    //ROS_WARN("[Bezier Trajectory] Start stacking the objective");
    
    int min_order_l = floor(minimize_order);
    int min_order_u = ceil (minimize_order);

    
    hessian.resize(numberOfVariables,numberOfVariables);

   
    {    
        int sub_shift = 0;
        
        for(int k = 0; k < segment_num; k ++)
        {
            double scale_k = corridor[k].t; //公式（7）的缩放因子 scale_k 在此处即为 每段轨迹所对应的时间 corridor[k].t
            for(int p = 0; p < 3; p ++ ) // x,y,z 三个方向
                for( int i = 0; i < s1d1CtrlP_num; i ++ )  //一段轨迹在一个方向上控制点数
                    for( int j = 0; j < s1d1CtrlP_num; j ++ )
                    {
                        int qsubi,qsubj;
                        double qval;
                        qsubi = sub_shift + p * s1d1CtrlP_num + i;   //qsubi[idx]对应 Q矩阵元素 行坐标
                        qsubj = sub_shift + p * s1d1CtrlP_num + j;   //qsubj[idx]对应 Q矩阵元素 列坐标
                        // qsubi = k*s1CtrlP_num + p * s1d1CtrlP_num + i;   //qsubi[idx]对应 Q矩阵元素 行坐标
                        // qsubj = k*s1CtrlP_num + p * s1d1CtrlP_num + j;   //qsubj[idx]对应 Q矩阵元素 列坐标

                        if(min_order_l == min_order_u)
                            qval  = MQM(i, j) /(double)pow(scale_k, 2 * min_order_u - 3); // qval[idx] 为对应 （qsubi[idx]，qsubj[idx] ）位置处的值
                        else                                                                   // pow(scale_k, 2 * min_order_u - 3) 系数对应 论文 公式（7）
                            qval = ( (minimize_order - min_order_l) / (double)pow(scale_k, 2 * min_order_u - 3)
                                        + (min_order_u - minimize_order) / (double)pow(scale_k, 2 * min_order_l - 3) ) * MQM(i, j); //加权平均

                        hessian.insert(qsubi,qsubj)=qval;
                    }
                       
            sub_shift += s1CtrlP_num; //s1CtrlP_num 一段轨迹 x y z 三个方向所对应控制点总数
        }
    }
         
    ros::Time time_end1 = ros::Time::now();

    Eigen::VectorXd QPSolution;

    // instantiate the solver
    OsqpEigen::Solver solver;

    // settings
    //solver.settings()->setVerbosity(false);
    solver.settings()->setWarmStart(true);

    // set the initial data of the QP solver
    solver.data()->setNumberOfVariables(numberOfVariables);
    solver.data()->setNumberOfConstraints(numberOfConstrains);
    if(!solver.data()->setHessianMatrix(hessian)) return -1;
    if(!solver.data()->setGradient(gradient)) return -1;
    if(!solver.data()->setLinearConstraintsMatrix(linearMatrix)) return -1;
    if(!solver.data()->setLowerBound(lowerBound)) return -1;
    if(!solver.data()->setUpperBound(upperBound)) return -1;

    // instantiate the solver
    if(!solver.initSolver()) return -1;

    //solve the QP problem
    // if(!solver.solve()) return 1;
    
    if(!solver.solve()){
      ROS_WARN("In solver, falied ");
      return -1;
    }
    // get the solution
    QPSolution = solver.getSolution();
    //std::cout<<QPSolution<<std::endl;
     
     

    ros::Time time_end2 = ros::Time::now();
    ROS_WARN("time consume in optimize is :");
    cout<<time_end2 - time_end1<<endl;

    

    VectorXd d_var(ctrlP_num);
    for(int i = 0; i < ctrlP_num; i++)
        d_var(i) = QPSolution[i];
    
    PolyCoeff = MatrixXd::Zero(segment_num, 3 *(traj_order + 1) ); //轨迹系数存储方式 每一行存储一条轨迹参数 顺序为 c_x c_y c_z ,共 segment_num 行

    int var_shift = 0;
    for(int i = 0; i < segment_num; i++ ) // 对应 PolyCoeff 的行
    {
        for(int j = 0; j < 3 * n_poly; j++) // 对应 PolyCoeff 的列  n_poly 为一段轨迹在一个方向上参数个数
            PolyCoeff(i , j) = d_var(j + var_shift);

        var_shift += 3 * n_poly;
    }   

    return 1;
}











// static void MSKAPI printstr(void *handle, MSKCONST char str[])
// {
//   printf("%s",str);
// }

// int TrajectoryGenerator::BezierPloyCoeffGeneration(
//             const vector<Cube> &corridor,
//             const MatrixXd &MQM,
//             const MatrixXd &pos,
//             const MatrixXd &vel,
//             const MatrixXd &acc,
//             const double maxVel,
//             const double maxAcc,
//             const int traj_order,
//             const double minimize_order,
//             const double margin,
//             const bool & isLimitVel,
//             const bool & isLimitAcc,
//             double & obj,
//             MatrixXd & PolyCoeff)  // define the order to which we minimize.   1 -- velocity, 2 -- acceleration, 3 -- jerk, 4 -- snap  
// {   
// #define ENFORCE_VEL  isLimitVel // whether or not adding extra constraints for ensuring the velocity feasibility
// #define ENFORCE_ACC  isLimitAcc // whether or not adding extra constraints for ensuring the acceleration feasibility

//     double initScale = corridor.front().t; //获取第一段轨迹时间
//     double lstScale  = corridor.back().t;
//     int segment_num  = corridor.size();

//     int n_poly = traj_order + 1;
//     int s1d1CtrlP_num = n_poly;  //一段轨迹在一个方向的控制点个数
//     int s1CtrlP_num   = 3 * s1d1CtrlP_num; //一段轨迹的总控制点个数

//     int equ_con_s_num = 3 * 3; // p, v, a in x, y, z axis at the start point
//     int equ_con_e_num = 3 * 3; // p, v, a in x, y, z axis at the end point
//     int equ_con_continuity_num = 3 * 3 * (segment_num - 1);  //中间点 p v a 连续性约束
//     int equ_con_num   = equ_con_s_num + equ_con_e_num + equ_con_continuity_num; // p, v, a in x, y, z axis in each segment's joint position
    
//     int vel_con_num = 3 *  traj_order * segment_num; //速度 控制点个数
//     int acc_con_num = 3 * (traj_order - 1) * segment_num;  //加速度 控制点个数

//     if( !ENFORCE_VEL )  
//         vel_con_num = 0;

//     if( !ENFORCE_ACC )
//         acc_con_num = 0;

//     int high_order_con_num = vel_con_num + acc_con_num; 
//     //int high_order_con_num = 0; //3 * traj_order * segment_num;

//     int con_num   = equ_con_num + high_order_con_num;  //总约束个数
//     int ctrlP_num = segment_num * s1CtrlP_num;  //所有段轨迹的控制点总数

//     double x_var[ctrlP_num]; //存储 MOSEK 最终优化结果
//     double primalobj;
    
//     /* con_bdk 存储每个约束的类型及约束值 */
//     MSKrescodee  r; 
//     vector< pair<MSKboundkeye, pair<double, double> > > con_bdk;

//     /*         Bound keys as defined in the enum MSKboundkeye
//       Bound key       Type of bound       Lower bound        Upper bound

//       MSK_BK_FX        u_j=l_j            finite             finite
//       MSK_BK_FR        free               -inf                +inf
//       MSK_BK_LO        l_j<=···           finite              +inf
//       MSK_BK_RA        l_j<=···<=u_j      finite             finite
//       MSK_BK_UP        ···<=u_j           -inf               finite

//     */
    
//     if(ENFORCE_VEL)
//     {
//         /***  Stack the bounding value for the linear inequality for the velocity constraints  ***/
//         for(int i = 0; i < vel_con_num; i++)
//         {
//             pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair( MSK_BK_RA, make_pair( - maxVel,  + maxVel) );
//             con_bdk.push_back(cb_ie);   
//         }
//     }

//     if(ENFORCE_ACC)
//     {
//         /***  Stack the bounding value for the linear inequality for the acceleration constraints  ***/
//         for(int i = 0; i < acc_con_num; i++)
//         {
//             pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair( MSK_BK_RA, make_pair( - maxAcc,  maxAcc) ); 
//             con_bdk.push_back(cb_ie);   
//         }
//     }

//     /*等式约束的顺序
//     [s_p_x, s_p_y, s_p_z, s_v_x, s_v_y, s_v_z, s_a_x, s_a_y, s_a_z, e_p_x, e_p_y, e_p_z, e_v_x, e_v_y, e_v_z, e_a_x, e_a_y, e_a_z, 0, 0, ··· ··· 0, 0}'
//      起点 位置 速度 加速度 等式约束
//      终点 位置 速度 加速度 等式约束
//      中间点 位置 速度 加速度 等式约束
//     */

//     //ROS_WARN("[Bezier Trajectory] equality bound %d", equ_con_num);
//     for(int i = 0; i < equ_con_num; i ++ ){ 
//         double beq_i;
//         if(i < 3)                    beq_i = pos(0, i);      //start point: x y z
//         else if (i >= 3  && i < 6  ) beq_i = vel(0, i - 3);  //start point: v_x v_y v_z
//         else if (i >= 6  && i < 9  ) beq_i = acc(0, i - 6);  //start point: a_x a_y a_z
//         else if (i >= 9  && i < 12 ) beq_i = pos(1, i - 9 ); //end point: x y z
//         else if (i >= 12 && i < 15 ) beq_i = vel(1, i - 12); //end point: v_x v_y v_z
//         else if (i >= 15 && i < 18 ) beq_i = acc(1, i - 15); //end point: a_x a_y a_z
//         else beq_i = 0.0;

//         pair<MSKboundkeye, pair<double, double> > cb_eq = make_pair( MSK_BK_FX, make_pair( beq_i, beq_i ) ); // # cb_eq means: constriants boundary of equality constrain
//         con_bdk.push_back(cb_eq); 
//     }

//     /* ## define a container for control points' boundary and boundkey ## */ 
//     /* ## dataType in one tuple is : boundary type, lower bound, upper bound ## */
//     vector< pair<MSKboundkeye, pair<double, double> > > var_bdk; 

//     for(int k = 0; k < segment_num; k++)  // segment_num 段轨迹
//     {   
//         Cube cube_     = corridor[k];  //获取每段轨迹所对应的走廊参数
//         double scale_k = cube_.t;    //scale factor 将每段轨迹时间映射到 [ 0,1] 区间

//         for(int i = 0; i < 3; i++ )  //三个方向 x y z 
//         {   
//             for(int j = 0; j < n_poly; j ++ ) // n_poly 个控制点
//             {   
//                 pair<MSKboundkeye, pair<double, double> > vb_x;

//                 double lo_bound, up_bound;
                
//                 if(k > 0)
//                 {
//                     lo_bound = (cube_.box[i].first  + margin) / scale_k;  // margin 立方体边距 
//                     up_bound = (cube_.box[i].second - margin) / scale_k;  //对应论文 公式（5） scale_k 本应在等式左侧，此处在等式两端同时除以 scale_k
//                 }
//                 else //第一个cube是基准
//                 {
//                     lo_bound = (cube_.box[i].first)  / scale_k;
//                     up_bound = (cube_.box[i].second) / scale_k;
//                 }

//                 vb_x  = make_pair( MSK_BK_RA, make_pair( lo_bound, up_bound ) ); // # vb_x means: varialbles boundary of unknowns x (Polynomial coeff)

//                 var_bdk.push_back(vb_x);
//             }
//         } 
//     }

//     MSKint32t  j,i; 
//     MSKenv_t   env; 
//     MSKtask_t  task; 

//     /***********************************************************MOSEK 轨迹优化求解部分*******************************************************************/
//     // Create the mosek environment. 
//     r = MSK_makeenv( &env, NULL ); 
  
//     // Create the optimization task. 
//     r = MSK_maketask(env,con_num, ctrlP_num, &task);  //参数：MSKenv 等式约束数目 待优化参数总数 MSKtask

// // Parameters used in the optimizer
// //######################################################################
//     //MSK_putintparam (task, MSK_IPAR_OPTIMIZER , MSK_OPTIMIZER_INTPNT );
//     MSK_putintparam (task, MSK_IPAR_NUM_THREADS, 1);
//     MSK_putdouparam (task, MSK_DPAR_CHECK_CONVEXITY_REL_TOL, 1e-2);
//     MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_DFEAS,  1e-4);
//     MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_PFEAS,  1e-4);
//     MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_INFEAS, 1e-4);
//     //MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_REL_GAP, 5e-2 );
// //######################################################################
    
//     //r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr); 
    
//     // Append empty constraints. 
//      //The constraints will initially have no bounds. 
//     if ( r == MSK_RES_OK ) 
//       r = MSK_appendcons(task,con_num);  

//     // Append optimizing variables. The variables will initially be fixed at zero (x=0). 
//     if ( r == MSK_RES_OK ) 
//       r = MSK_appendvars(task,ctrlP_num); 
    
//     //ROS_WARN("set variables boundary");
//     for(j = 0; j<ctrlP_num && r == MSK_RES_OK; ++j){ 
//         if (r == MSK_RES_OK) 
//             r = MSK_putvarbound(task, 
//                                 j,                            // Index of variable. 
//                                 var_bdk[j].first,             // Bound key.
//                                 var_bdk[j].second.first,      // Numerical value of lower bound.
//                                 var_bdk[j].second.second );   // Numerical value of upper bound.      
//     } 
    
//     // Set the bounds on constraints.  此处对应 论文 公式（11）
//     //   for i=1, ...,con_num : blc[i] <= constraint i <= buc[i] 
//     for( i = 0; i < con_num && r == MSK_RES_OK; i++ ) {
//         r = MSK_putconbound(task, 
//                             i,                            // Index of constraint. 
//                             con_bdk[i].first,             // Bound key.
//                             con_bdk[i].second.first,      // Numerical value of lower bound.
//                             con_bdk[i].second.second );   // Numerical value of upper bound. 
//     }

//     /***************************  Matrix A *************************************/
//     /**************************   AC<=B  ***************************************/
    
//     //ROS_WARN("[Bezier Trajectory] Start stacking the Linear Matrix A, inequality part");
//     int row_idx = 0; //Matrix A 的行号
//     // The velocity constraints
//     if(ENFORCE_VEL)
//     {   
//         for(int k = 0; k < segment_num ; k ++ )
//         {   
//             for(int i = 0; i < 3; i++)
//             {  // for x, y, z loop
//                 for(int p = 0; p < traj_order; p++)
//                 {
//                     int nzi = 2;
//                     MSKint32t asub[nzi];
//                     double aval[nzi];

//                     aval[0] = -1.0 * traj_order;  //此处对应 论文 公式（12）
//                     aval[1] =  1.0 * traj_order;

//                     asub[0] = k * s1CtrlP_num + i * s1d1CtrlP_num + p;    // 此处对应约束在MOSEK中的序号
//                     asub[1] = k * s1CtrlP_num + i * s1d1CtrlP_num + p + 1;    //s1CtrlP_num :一段轨迹的控制点数 = 3* s1d1CtrlP_num：一段轨迹在一个方向的控制点数 
  
//                     r = MSK_putarow(task, row_idx, nzi, asub, aval);   //将 Matrix A 的 第 row_idx 行 第 asub 起的 nzi 个数 替换为 aval 中的数
//                     row_idx ++;
//                 }
//             }
//         }
//     }

//     // The acceleration constraints
//     if(ENFORCE_ACC)
//     {
//         for(int k = 0; k < segment_num ; k ++ )
//         {
//             for(int i = 0; i < 3; i++)
//             { 
//                 for(int p = 0; p < traj_order - 1; p++)
//                 {    
//                     int nzi = 3;
//                     MSKint32t asub[nzi];
//                     double aval[nzi];

//                     aval[0] =  1.0 * traj_order * (traj_order - 1) / corridor[k].t;  //此处对应 论文 公式（12）
//                     aval[1] = -2.0 * traj_order * (traj_order - 1) / corridor[k].t;
//                     aval[2] =  1.0 * traj_order * (traj_order - 1) / corridor[k].t;
//                     asub[0] = k * s1CtrlP_num + i * s1d1CtrlP_num + p;    
//                     asub[1] = k * s1CtrlP_num + i * s1d1CtrlP_num + p + 1;    
//                     asub[2] = k * s1CtrlP_num + i * s1d1CtrlP_num + p + 2;    
                    
//                     r = MSK_putarow(task, row_idx, nzi, asub, aval);    //将 Matrix A 的 第 row_idx 行 第 asub 起的 nzi 个数 替换为 aval 中的数
//                     row_idx ++;
//                 }
//             }
//         }
//     }
//     /*   Start position  */
//     {
//         // position :
//         for(int i = 0; i < 3; i++)
//         {  // loop for x, y, z       
//             int nzi = 1;
//             MSKint32t asub[nzi];
//             double aval[nzi];
//             aval[0] = 1.0 * initScale;
//             asub[0] = i * s1d1CtrlP_num;
//             r = MSK_putarow(task, row_idx, nzi, asub, aval);    
//             row_idx ++;
//         }
//         // velocity :
//         for(int i = 0; i < 3; i++)
//         {  // loop for x, y, z       
//             int nzi = 2;
//             MSKint32t asub[nzi];
//             double aval[nzi];
//             aval[0] = - 1.0 * traj_order;
//             aval[1] =   1.0 * traj_order;
//             asub[0] = i * s1d1CtrlP_num;
//             asub[1] = i * s1d1CtrlP_num + 1;
//             r = MSK_putarow(task, row_idx, nzi, asub, aval);   
//             row_idx ++;
//         }
//         // acceleration : 
//         for(int i = 0; i < 3; i++)
//         {  // loop for x, y, z       
//             int nzi = 3;
//             MSKint32t asub[nzi];
//             double aval[nzi];
//             aval[0] =   1.0 * traj_order * (traj_order - 1) / initScale;
//             aval[1] = - 2.0 * traj_order * (traj_order - 1) / initScale;
//             aval[2] =   1.0 * traj_order * (traj_order - 1) / initScale;
//             asub[0] = i * s1d1CtrlP_num;
//             asub[1] = i * s1d1CtrlP_num + 1;
//             asub[2] = i * s1d1CtrlP_num + 2;
//             r = MSK_putarow(task, row_idx, nzi, asub, aval);    
//             row_idx ++;
//         }
//     }      

//     /*   End position  */
//     //ROS_WARN(" end position");
//     {   
//         // position :
//         for(int i = 0; i < 3; i++)
//         {  // loop for x, y, z       
//             int nzi = 1;
//             MSKint32t asub[nzi];
//             double aval[nzi];
//             asub[0] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num;
//             aval[0] = 1.0 * lstScale;
//             r = MSK_putarow(task, row_idx, nzi, asub, aval);    
//             row_idx ++;
//         }
//         // velocity :
//         for(int i = 0; i < 3; i++)
//         { 
//             int nzi = 2;
//             MSKint32t asub[nzi];
//             double aval[nzi];
//             asub[0] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num - 1;
//             asub[1] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num;
//             aval[0] = - 1.0;
//             aval[1] =   1.0;
//             r = MSK_putarow(task, row_idx, nzi, asub, aval);    
//             row_idx ++;
//         }
//         // acceleration : 
//         for(int i = 0; i < 3; i++)
//         { 
//             int nzi = 3;
//             MSKint32t asub[nzi];
//             double aval[nzi];
//             asub[0] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num - 2;
//             asub[1] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num - 1;
//             asub[2] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num;
//             aval[0] =   1.0 / lstScale;
//             aval[1] = - 2.0 / lstScale;
//             aval[2] =   1.0 / lstScale;
//             r = MSK_putarow(task, row_idx, nzi, asub, aval);    
//             row_idx ++;
//         }
//     }

//     /*   joint points  */
//     /********************************************** 对应论文 公式（10）**********************************************/
//     //ROS_WARN(" joint position");
//     {
//         int sub_shift = 0;
//         double val0, val1;
//         for(int k = 0; k < (segment_num - 1); k ++ )
//         {   
//             double scale_k = corridor[k].t;
//             double scale_n = corridor[k+1].t;
//             // position :
//             val0 = scale_k;
//             val1 = scale_n;
//             for(int i = 0; i < 3; i++)  
//             {  // loop for x, y, z
//                 int nzi = 2;
//                 MSKint32t asub[nzi];
//                 double aval[nzi];

//                 // This segment's last control point
//                 aval[0] = 1.0 * val0;
//                 asub[0] = sub_shift + (i+1) * s1d1CtrlP_num - 1;

//                 // Next segment's first control point
//                 aval[1] = -1.0 * val1;
//                 asub[1] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num;
//                 r = MSK_putarow(task, row_idx, nzi, asub, aval);    
//                 row_idx ++;
//             }
//             // velocity
//             // val0 = pow(scale_k,1-1)=1
//             // val1 = pow(scale_n,1-1)=1
//             for(int i = 0; i < 3; i++)
//             {  
//                 int nzi = 4;
//                 MSKint32t asub[nzi];
//                 double aval[nzi];
                
//                 // This segment's last velocity control point
//                 aval[0] = -1.0;  
//                 aval[1] =  1.0;
//                 asub[0] = sub_shift + (i+1) * s1d1CtrlP_num - 2;    
//                 asub[1] = sub_shift + (i+1) * s1d1CtrlP_num - 1;   
//                 // Next segment's first velocity control point
//                 aval[2] =  1.0;  
//                 aval[3] = -1.0;

//                 asub[2] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num;    
//                 asub[3] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num + 1;

//                 r = MSK_putarow(task, row_idx, nzi, asub, aval);    
//                 row_idx ++;
//             }
//             // acceleration :
//             val0 = 1.0 / scale_k;
//             val1 = 1.0 / scale_n;
//             for(int i = 0; i < 3; i++)
//             {  
//                 int nzi = 6;
//                 MSKint32t asub[nzi];
//                 double aval[nzi];
                
//                 // This segment's last velocity control point
//                 aval[0] =  1.0  * val0;
//                 aval[1] = -2.0  * val0;
//                 aval[2] =  1.0  * val0;
//                 asub[0] = sub_shift + (i+1) * s1d1CtrlP_num - 3;    
//                 asub[1] = sub_shift + (i+1) * s1d1CtrlP_num - 2;   
//                 asub[2] = sub_shift + (i+1) * s1d1CtrlP_num - 1;   
//                 // Next segment's first velocity control point
//                 aval[3] =  -1.0  * val1;
//                 aval[4] =   2.0  * val1;
//                 aval[5] =  -1.0  * val1;
//                 asub[3] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num;    
//                 asub[4] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num + 1;
//                 asub[5] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num + 2;

//                 r = MSK_putarow(task, row_idx, nzi, asub, aval);    
//                 row_idx ++;
//             }

//             sub_shift += s1CtrlP_num;
//         }
//     }

//     //ROS_WARN("[Bezier Trajectory] Start stacking the objective");
    
//     int min_order_l = floor(minimize_order);
//     int min_order_u = ceil (minimize_order);

//     int NUMQNZ = 0;
//     for(int i = 0; i < segment_num; i ++)
//     {   
//         int NUMQ_blk = (traj_order + 1);                       // default minimize the jerk and minimize_order = 3
//         NUMQNZ      += 3 * NUMQ_blk * (NUMQ_blk + 1) / 2;      //计算Q矩阵下三角中对应元素的个数 Q：（traj_order + 1）*（traj_order + 1）维 
//     }                                                          // 1+2+3+······+（traj_order + 1）=（1+traj_order + 1）*（traj_order + 1）/2
//     MSKint32t  qsubi[NUMQNZ], qsubj[NUMQNZ];
//     double     qval[NUMQNZ];
    
//     /*
//     The lower triangular part of the matrix Q is specified using an unordered sparse triplet format（无序稀疏三元组格式）
//     1： only non-zero elements are specified (any element not specified is 0 by definition),
//     2： the order of the non-zero elements is insignificant
//     3： only the lower triangular part should be specified
//     */
//     {    
//         int sub_shift = 0;
//         int idx = 0;
//         for(int k = 0; k < segment_num; k ++)
//         {
//             double scale_k = corridor[k].t; //公式（7）的缩放因子 scale_k 在此处即为 每段轨迹所对应的时间 corridor[k].t
//             for(int p = 0; p < 3; p ++ ) // x,y,z 三个方向
//                 for( int i = 0; i < s1d1CtrlP_num; i ++ )  //一段轨迹在一个方向上控制点数
//                     for( int j = 0; j < s1d1CtrlP_num; j ++ )
//                         if( i >= j )  // min : C‘QC  矩阵Q为对称矩阵（若Q为非对称矩阵，可通过变换转化为对称阵），因此只给 Q的下三角部分赋值
//                         {
//                             qsubi[idx] = sub_shift + p * s1d1CtrlP_num + i;   //qsubi[idx]对应 Q矩阵元素 行坐标
//                             qsubj[idx] = sub_shift + p * s1d1CtrlP_num + j;   //qsubj[idx]对应 Q矩阵元素 列坐标
//                             //qval[idx]  = MQM(i, j) /(double)pow(scale_k, 3);
//                             if(min_order_l == min_order_u)
//                                 qval[idx]  = MQM(i, j) /(double)pow(scale_k, 2 * min_order_u - 3); // qval[idx] 为对应 （qsubi[idx]，qsubj[idx] ）位置处的值
//                             else                                                                   // pow(scale_k, 2 * min_order_u - 3) 系数对应 论文 公式（7）
//                                 qval[idx] = ( (minimize_order - min_order_l) / (double)pow(scale_k, 2 * min_order_u - 3)
//                                             + (min_order_u - minimize_order) / (double)pow(scale_k, 2 * min_order_l - 3) ) * MQM(i, j); //加权平均
//                             idx ++ ;
//                         }

//             sub_shift += s1CtrlP_num; //s1CtrlP_num 一段轨迹 x y z 三个方向所对应控制点总数
//         }
//     }
         
//     ros::Time time_end1 = ros::Time::now();

//     //The quadratic objective is specified using the function MSK_putqobj（）
//     if ( r== MSK_RES_OK )
//          r = MSK_putqobj(task,NUMQNZ,qsubi,qsubj,qval);  //Replaces all quadratic terms in the objective.

//     /*MINIMIZE objective function. */
//     if ( r==MSK_RES_OK ) 
//          r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE); //Sets the objective sense.
    
//     //ros::Time time_opt = ros::Time::now();
//     bool solve_ok = false;
//     if ( r==MSK_RES_OK ) 
//       { 
//         //ROS_WARN("Prepare to solve the problem ");   
//         MSKrescodee trmcode; 

//         /* Run optimizer */
//         r = MSK_optimizetrm(task,&trmcode); 

//         /* Print a summary containing information
//         about the solution for debugging purposes. */
//         MSK_solutionsummary (task,MSK_STREAM_LOG); 
          
//         if ( r==MSK_RES_OK ) 
//         { 
//           /*Extracting the solution.
//           After optimizing the status of the solution is examined with a call to MSK_getsolsta（）
//           */
//           MSKsolstae solsta; 
//           MSK_getsolsta (task,MSK_SOL_ITR,&solsta); 
           
//           switch(solsta) 
//           { 
//             case MSK_SOL_STA_OPTIMAL:    
//             case MSK_SOL_STA_NEAR_OPTIMAL: 
              
//             /*The MSK_getxx function obtains the solution*/

//             r = MSK_getxx(task, 
//                           MSK_SOL_ITR,    // Request the interior solution.  
//                           x_var); 
            
//             /*After the meaning and quality of the solution (or certificate) have been established, 
//             we can query for the actual numerical values.*/
//             r = MSK_getprimalobj(  //the primal and dual objective value
//                 task,
//                 MSK_SOL_ITR,
//                 &primalobj);

//             obj = primalobj;
//             solve_ok = true;
            
//             break; 
            
//             case MSK_SOL_STA_DUAL_INFEAS_CER: 
//             case MSK_SOL_STA_PRIM_INFEAS_CER: 
//             case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER: 
//             case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:   
//               printf("Primal or dual infeasibility certificate found.\n"); 
//               break; 
               
//             case MSK_SOL_STA_UNKNOWN: 
//               printf("The status of the solution could not be determined.\n"); 
//               //solve_ok = true; // debug
//               break; 
//             default: 
//               printf("Other solution status."); 
//               break; 
//           } 
//         } 
//         else 
//         { 
//           printf("Error while optimizing.\n"); 
//         } 
//       }
     
//       if (r != MSK_RES_OK) 
//       { 
//         // In case of an error print error code and description. 
//         char symname[MSK_MAX_STR_LEN]; 
//         char desc[MSK_MAX_STR_LEN]; 
         
//         printf("An error occurred while optimizing.\n");      
//         MSK_getcodedesc (r, 
//                          symname, 
//                          desc); 
//         printf("Error %s - '%s'\n",symname,desc); 
//       } 
    
//     /* Delete the task and the associated data. */
//     MSK_deletetask(&task);

//     /* Delete the environment and the associated data. */ 
//     MSK_deleteenv(&env); 

//     ros::Time time_end2 = ros::Time::now();
//     ROS_WARN("time consume in optimize is :");
//     cout<<time_end2 - time_end1<<endl;

//     if(!solve_ok){
//       ROS_WARN("In solver, falied ");
//       return -1;
//     }

//     VectorXd d_var(ctrlP_num);
//     for(int i = 0; i < ctrlP_num; i++)
//         d_var(i) = x_var[i];
    
//     PolyCoeff = MatrixXd::Zero(segment_num, 3 *(traj_order + 1) ); //轨迹系数存储方式 每一行存储一条轨迹参数 顺序为 c_x c_y c_z ,共 segment_num 行

//     int var_shift = 0;
//     for(int i = 0; i < segment_num; i++ ) // 对应 PolyCoeff 的行
//     {
//         for(int j = 0; j < 3 * n_poly; j++) // 对应 PolyCoeff 的列  n_poly 为一段轨迹在一个方向上参数个数
//             PolyCoeff(i , j) = d_var(j + var_shift);

//         var_shift += 3 * n_poly;
//     }   

//     return 1;
// }
