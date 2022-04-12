#ifndef _DATA_TYPE_
#define _DATA_TYPE_

#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>
#include <eigen3/Eigen/Dense>
#include <vector>

#define inf 1>>30

struct Cube;
struct GridNode;
typedef GridNode* GridNodePtr;

struct Cube
{     
      //Eigen::Vector3d p1, p2, p3, p4, p5, p6, p7, p8;   // the 8 vertex of a cube 
      Eigen::MatrixXd vertex;   //8*3矩阵
      Eigen::Vector3d center; // the center of the cube
      bool valid;    // indicates whether this cube should be deleted

      double t; // time allocated to this cube
      std::vector< std::pair<double, double> > box;
/*
           P4------------P3 
           /|           /|              ^
          / |          / |              | z
        P1--|---------P2 |              |
         |  P8--------|--p7             |
         | /          | /               /--------> y
         |/           |/               /  
        P5------------P6              / x
*/                                                                                 

      // create a cube using 8 vertex and the center point
      Cube( Eigen::MatrixXd vertex_, Eigen::Vector3d center_)
      {
            vertex = vertex_;
            center = center_;
            valid = true;
            t = 0.0;
            box.resize(3);
      }

      // create a inscribe cube of a ball using the center point and the radius of the ball 内切于球
      void setVertex( Eigen::MatrixXd vertex_, double resolution_)
      {     
            vertex = vertex_;
            //p1 p4 p8 p5 左侧平面
            vertex(0,1) -= resolution_ / 2.0;  //p1_y
            vertex(3,1) -= resolution_ / 2.0;  //p4_y
            vertex(4,1) -= resolution_ / 2.0;  //p5_y
            vertex(7,1) -= resolution_ / 2.0;  //p8_y
            
            //p2 p3 p7 p6 右侧平面
            vertex(1,1) += resolution_ / 2.0;  //p2_y
            vertex(2,1) += resolution_ / 2.0;  //p3_y
            vertex(5,1) += resolution_ / 2.0;  //p6_y
            vertex(6,1) += resolution_ / 2.0;  //p7_y
            
            //p1 p2 p5 p6 正前平面
            vertex(0,0) += resolution_ / 2.0;  //p1_x
            vertex(1,0) += resolution_ / 2.0;  //p2_x
            vertex(4,0) += resolution_ / 2.0;  //p5_x
            vertex(5,0) += resolution_ / 2.0;  //p6_x
            
            //p4 p3 p7 p8 正后平面
            vertex(2,0) -= resolution_ / 2.0;  //p3_x
            vertex(3,0) -= resolution_ / 2.0;  //p4_x
            vertex(6,0) -= resolution_ / 2.0;  //p7_x
            vertex(7,0) -= resolution_ / 2.0;  //p8_x
            
            //p1 p2 p3 p4 正上平面
            vertex(0,2) += resolution_ / 2.0;  //p1_z
            vertex(1,2) += resolution_ / 2.0;  //p2_z
            vertex(2,2) += resolution_ / 2.0;  //p3_z
            vertex(3,2) += resolution_ / 2.0;  //p4_z
            
            //p5 p6 p7 p8 正下平面
            vertex(4,2) -= resolution_ / 2.0;  //p5_z
            vertex(5,2) -= resolution_ / 2.0;  //p6_z
            vertex(6,2) -= resolution_ / 2.0;  //p7_z
            vertex(7,2) -= resolution_ / 2.0;  //p8_z
            
            setBox();
      }
      
      void setBox()
      {
            box.clear();
            box.resize(3);
            box[0] = std::make_pair( vertex(3, 0), vertex(0, 0) );  //x轴向范围  顺序：小 大
            box[1] = std::make_pair( vertex(0, 1), vertex(1, 1) );  //y轴向范围
            box[2] = std::make_pair( vertex(4, 2), vertex(1, 2) );  //z轴向范围
      }

      void printBox()
      {
            std::cout<<"center of the cube: \n"<<center<<std::endl;
            std::cout<<"vertex of the cube: \n"<<vertex<<std::endl;
      }

      Cube()
      {  
         center = Eigen::VectorXd::Zero(3);
         vertex = Eigen::MatrixXd::Zero(8, 3);

         valid = true;
         t = 0.0;
         box.resize(3);
      }

      ~Cube(){}
};

struct GridNode
{     
   int id;        // 1--> open set, -1 --> closed set
   Eigen::Vector3d coord;
   Eigen::Vector3i index;
   
   double gScore, fScore;
   GridNodePtr cameFrom;
   std::multimap<double, GridNodePtr>::iterator nodeMapIt;
   double occupancy; 

   std::vector<GridNodePtr> hisNodeList; // use a list to record nodes in its history

   GridNode(Eigen::Vector3i _index)
   {  
      id = 0;
      index = _index;
      
      gScore = inf;
      fScore = inf;
      cameFrom = NULL;
   }

   GridNode(Eigen::Vector3i _index, Eigen::Vector3d _coord)
   {  
      id = 0;
      index = _index;
      coord = _coord;

      gScore = inf;
      fScore = inf;
      cameFrom = NULL;
   }

   GridNode(){};
   
   ~GridNode(){};
};

#endif