//Lara Querciagrossa
#include <stdio.h>
#include <cstdlib>
#include <memory>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <strstream>
#include <cmath>
#include "math.hpp"
#include "math-quat.hpp"

#define THIRD 0.333333333333333
#define ROOTTHREE 1.73205080756888       
#define PI 3.141592653589793238462643383279502884197
 
using namespace std;


int main(int argc, char *argv[])
{

  if (argc != 2)
    {
      std::cout << "[ERROR] Wrong number of input parameters: ./oxDNA_expand <oxDNAdump> " << std::endl;
      return 0;
    }

  if (string(argv[1]) == "--help" || argv[1] == "-h") {
    std::cout << " " << std::endl;
    std::cout << " Welcome to the oxDNA_exapand code!" << std::endl;
    std::cout << " This code will convert an oxDNA dump file (obtained with LAMMPS) contained only center of nucleotide" << std::endl;
    std::cout << " into a dump file with one bead in the backbone area and another one representing the base itself" << std::endl;
    std::cout << " " << std::endl;
    std::cout << " ./oxDNA_expand <oxDNAdump>" << std::endl;
  }
  else {

    ifstream conf ( argv[1] );
    if ( !conf.is_open() )
      {
	cout<<"ERROR: Could not open dump file\n";
	return 0;
      }

    std::cout << "Your configuration file is: " << argv[1] << std::endl;
    std::cout << " " << std::endl;

   /** Read configuration(s) **/

    //Count lines in files
    string line; 
    int totlines = 0;
    int totsteps = 0;
    string itemtimestep = "ITEM: TIMESTEP";
    while ( getline(conf,line) )
      {
	string first = line;
	if ( first.compare(itemtimestep) == 0 )
	  {
	    totsteps++;
	  }
	totlines++;
      }
    //conf.close();
    //std::cout << "Total lines in your configuration file: " << totlines << std::endl;
    std::cout << "In your configuration file, there are " << totsteps << " steps." << std::endl;

    //Configuration file has been read till the end: this position is cleared and than we set the new position to be read next at the beginning of the file
    conf.clear();
    conf.seekg(0, std::ios::beg);

    //Write an expanded dump file to be read in ovito
    ofstream fileout("oxDNAexpanded.dump", ios::out);

    string null;

    //File will be read one line at a time to save memory
    for (int i=0; i<totsteps; i++)
      //for (int i=0; i<1; i++)
      {
	std::cout << "Converting step " << i+1 << "/" << totsteps << std::endl;

	int timestep, atoms;
	double xlo, xhi, ylo, yhi, zlo, zhi;

	conf >> null >> null;
	conf >> timestep;
	conf >> null >> null >>  null >> null;
	conf >> atoms;
	conf >> null >> null >>  null >> null >> null >> null;
	conf >> xlo >> xhi;
	conf >> ylo >> yhi;
	conf >> zlo >> zhi;


	//Saving "start" part of dump file
	fileout << "ITEM: TIMESTEP" << std::endl;
	fileout << timestep << std::endl;
	fileout << "ITEM: NUMBER OF ATOMS" << std::endl;
	fileout << 3*atoms << std::endl;
	fileout << "ITEM: BOX BOUNDS" << std::endl;
	fileout << xlo << " " << xhi << std::endl;
	fileout << ylo << " " << yhi << std::endl;
	fileout << zlo << " " << zhi << std::endl;



	conf >> null >> null >> null >> null >> null >> null >> null >> null >> null >> null >> null >> null >> null >> null >> null >> null >> null >> null >> null >> null >> null;

	int cnt = 1;
	//Saving positions into dump file
	fileout << "ITEM: ATOMS id type x y z c_q[1] c_q[2] c_q[3] c_q[4] c_shape[1] c_shape[2] c_shape[3] vx vy vz angmomx angmomy angmomz mol" << std::endl;
	for (int b=0; b<atoms; b++)
	  {
	    int mol, id, type; //MoleculeID, AtomID, AtomType
	    double rx, ry, rz; //Position of com
	    double qw, qx, qy, qz; //Orientation
	    double sx, sy, sz; //Shape
	    double vx, vy, vz; //Velocities 
	    double wx, wy, wz; //Angular velocities

	    double bx, by, bz; //Backbone-baseversor
	    double nx, ny, nz; //Normalversor 

	    conf >> id >> type >> rx >> ry >> rz >> qw >> qx >> qy >> qz >> sx >> sy >> sz >> vx >> vy >> vz >> wx >> wy >> wz >> mol;

	    double quat[4] = {qw, qx, qy, qz};
	    double** mat = QuaternionToMatrix(quat);

	    bx = mat[0][0];
	    by = mat[0][1];
	    bz = mat[0][2];
	    nx = mat[2][0];
	    ny = mat[2][1];
	    nz = mat[2][2];

	    //Backbone site: (center) - 0.40 * (axis_vector)
	    fileout << cnt << " " << " 7 " << " " << rx-0.4*bx  << " " << ry-0.4*by  << " " << rz-0.4*bz <<  " 1.0 0.0 0.0 0.0 " << " 0.3 0.3 0.3 "  << " " << vx << " " << vy << " " << vz << " " << wx << " " << wy << " " << wz << " " << mol << std::endl; 
	    cnt++;

	    //Stacking site: (center) + 0.34 * (axis vector)
	    fileout << cnt << " " << " 9 " << " " << rx+0.34*bx  << " " << ry+0.34*by  << " " << rz+0.34*bz <<  " 1.0 0.0 0.0 0.0 " << " 0.1 0.1 0.1 "  << " " << vx << " " << vy << " " << vz << " " << wx << " " << wy << " " << wz << " " << mol << std::endl; 
	    cnt++;

	    //Base site: (center) + 0.40 * (axis vector)
	    fileout << cnt << " " << type << " " << rx+0.4*bx  << " " << ry+0.4*by  << " " << rz+0.4*bz <<  " " << qw << " " << qx << " " << qy << " " << qz << " " <<  " 0.4 0.3 0.1 " << " " << vx << " " << vy << " " << vz << " " << wx << " " << wy << " " << wz << " " << mol << std::endl; 
	    cnt++;
	  }

      }

    std::cout << " " << std::endl;
    std::cout << "Files oxDNAconv_expand.dump has been written." << std::endl;
    std::cout << " " << std::endl;

  }
}
