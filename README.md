# Btraj_osqp-eigen
Change the QP solver from Mosek to osqp-eigen ,make it possible to run it on you arm linux

# install [armadillo]:   
 sudo apt-get install libarmadillo-dev

# install OSQP :  
  git clone --recursive https://github.com/osqp/osqp  
  cd osqp  
  mkdir build  
  cd build  
  cmake -G "Unix Makefiles" ..  
  cmake --build .  
  cmake --build . --target install  
  
# install osqp-eigen  
  git clone https://github.com/robotology/osqp-eigen.git  
  cd osqp-eigen  
  mkdir build  
  cd build  
  cmake ../  
  make  
  make install  
  
# Clone the repository to your catkin workspace and catkin_make. For example:  
  cd ~/catkin_ws/src  
  git clone git@github.com:YYB666/Btraj_osqp-eigen.git  
  cd ../  
  catkin_make  
  source ~/catkin_ws/devel/setup.bash  
  
# Usage  
  roslaunch bezier_planer simulation.launch  
