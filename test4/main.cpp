#include <windows.h>
#include <gl/gl.h>
#include <gl/glu.h>
#include <GL/glui.h>
#include <stdio.h>
#include <GL/glut.h>
#include <math.h>
#include <vector>
#include <cmath>

using namespace std;

//struct to store the RBG value of a color
struct Color{
	double red;
	double blue;
	double green;
};

//struct that stores all the data pertaining to an individual particle
struct Particle{
	double x;
	double y;
	double z;

	double vx;
	double vy;
	double vz;

	Color color;

	double mass;
	double charge;
	double radius;

	double netForcex;
	double netForcey;
	double netForcez;

	double displacementx;
	double displacementy;
	double displacementz;
};
//pretty sure this is pie
const float PI = 4.0*atan(1.0);
double REALLY_SMALL = 0.000000000001;
int sphereQuality = 10;
double dt = 0.01;
double earthDisplacement = -1.0;
double G = 6.67 *pow(10.0, -11.0);
double ke = 8980000000;

double boxHeight = -.2;
double boxSize = 1.7;

float x,y,z = 0;
float worldX, worldY, worldZ = 0;
float centerX, centerY, centerZ = 0;
float angle = 0;
float rotation = 0.0;
float mouse1 = 0.0;
float mouse2 = 0.0;
float mouseX;
float mouseY;
float speed = 5.0;

bool earthOn = true;
bool animationOn = true;
bool collisionsOn = true;
bool boxOn = false;
bool newBall = true;
bool bBall = true;

//decided I wanted the "where" of the bBall hoops to be global vars
double xDistance = 3;
double boardSize = 3;
double boardHeight = 1.5;
double hoopRadius = 1;
double hoopThickness = .1;

//store all the particles in here
vector<Particle> particleWorld;

//use this to keep track of our main window
int mainWindow;



struct Paritcle particle(double mass, double charge, double radius, double x, double y, double z, double vx, double vy, double vz, int index);
void calcForces(double dt); //computes the forces between the particles
Particle collision(Particle currentParticle, Particle compParticle,double distanceBetween,int i,int j);

/*
*this is called every time OpenGL idles. Think of it as the "Next Frame" function.
*/
void Idle()
{

	if (animationOn){
		calcForces(dt);
	}

	glutSetWindow(mainWindow);
	//plan to call the calc forces stuff here
	glutPostRedisplay();
	//Sleep (100);
}

//constructor for a particle
struct Particle particle(double mass, double charge, double radius, double x, double y, double z, double vx, double vy, double vz, Color color){
	//when we make a point, we set it to location x,y,z with no acceleration, and velocity vx,vy,vz and mass mass, charge charge, radius radius
	Particle particle;


	particle.x = x; //x position
	particle.y = y; //y position
	particle.z = z; //z position
	particle.vx = vx;
	particle.vy = vy;
	particle.vz = vz;
	particle.color = color;
	particle.mass = mass;
	particle.charge = charge;
	particle.radius = radius;

	particleWorld.push_back(particle);

	return particle;
}


//finds the new locations for all the particles in the next time step
void calcForces(double dt)
{
	for (int i = 0 ; i < particleWorld.size() ; i++){
		Particle currentParticle = particleWorld[i];
		//reseting the net forces and displacement
		currentParticle.netForcex = 0;
		currentParticle.netForcey = 0;
		currentParticle.netForcez = 0;
		currentParticle.displacementx = 0;
		currentParticle.displacementy = 0;
		currentParticle.displacementz = 0;

		//adjust net Force if earth is on
		if (earthOn){
			currentParticle.netForcey = -9.8*currentParticle.mass;
			//the rest is looking redundant
			/*			//going to check for collision issues here
			if(currentParticle.y - currentParticle.radius <= earthDisplacement){
				//flip the verticle velocity
				currentParticle.vy = 0 - currentParticle.vy;
				//put the particle on the ground
				currentParticle.y = earthDisplacement + .00001 + currentParticle.radius;
				
			}
			*/
		}
		particleWorld[i] = currentParticle;
	}

	//Dealing with forces and collisions
	//calculate the force on i from j
	//going to allow massless point charges, so if mass of i = 0, we just set the force to 0
	for (int i = 0 ; i < particleWorld.size() ; i++){
		Particle currentParticle = particleWorld[i];

		for (int j = 0 ; j < particleWorld.size() ; j++){
			if(currentParticle.mass != 0 && i != j)
			{ // if it is massless, then we dont care, and there is no force if we are looking at the same particle
				Particle compParticle = particleWorld[j];

				double delta_x = pow(currentParticle.x - compParticle.x,2);
				double delta_y = pow(currentParticle.y - compParticle.y,2);
				double delta_z = pow(currentParticle.z - compParticle.z,2);
				double rSquared = delta_x + delta_y + delta_z;
				double distanceBetween = sqrt(rSquared);

				double Fg = G * currentParticle.mass * compParticle.mass / rSquared;
				double Fe = -ke * currentParticle.charge * compParticle.charge / rSquared;

				if(delta_x != 0){ //seem to randomly get 'nan' if we dont do this
					double Fx = (Fg + Fe) * sqrt(rSquared - delta_y - delta_z);
					if (currentParticle.x > compParticle.x) {
						Fx = Fx * (-1); //getting the direction right
					}
					currentParticle.netForcex += Fx;
				}

				if(delta_y !=0){
					double Fy = (Fg + Fe) * sqrt(rSquared - delta_x - delta_z);
					if (currentParticle.y > compParticle.y){
						Fy = Fy * (-1); //getting the direction right
					}
					currentParticle.netForcey += Fy;
				}

				if(delta_z != 0){
					double Fz = (Fg + Fe) * sqrt(rSquared - delta_y - delta_x);
					if (currentParticle.z > compParticle.z){
						Fz = Fz * (-1); //getting the direction right
					}
					currentParticle.netForcez += Fz;
				}
				
				if (collisionsOn){
					//decided to clean this up, all the math is taken care of elsewhere
					currentParticle = collision(currentParticle, compParticle,distanceBetween, i, j);
				}

				
			}
		}
				


		//lastly, we have to save this newly calculated data
		particleWorld[i] = currentParticle;

	}
				


	//now we move the particles accordingly
	for (int i = 0 ; i < particleWorld.size() ; i++){
		Particle currentParticle = particleWorld[i];
		if (currentParticle.mass != 0){


			//calculate delta v, change the speed accordingly
			//new velocity = old velocity + netforce * dt / mass
			currentParticle.vx += currentParticle.netForcex * dt / currentParticle.mass;
			currentParticle.vy += currentParticle.netForcey * dt / currentParticle.mass;
			currentParticle.vz += currentParticle.netForcez * dt / currentParticle.mass;


			//using the speed in theWorld, calculate the new position, and displacement
			//displacement is the old velocity * dt
			currentParticle.displacementx = currentParticle.vx*dt;
			currentParticle.displacementy = currentParticle.vy*dt;
			currentParticle.displacementz = currentParticle.vz*dt;


			//new position = old position + displacement
			currentParticle.x += currentParticle.displacementx;
			currentParticle.y += currentParticle.displacementy;
			currentParticle.z += currentParticle.displacementz;

			//if the earth is on, we see if the particle went through the ground this time step
			//if it did, we fix the velocity and position accordingly
			if (earthOn){
			if(currentParticle.y - currentParticle.radius <= earthDisplacement){
				//flip the verticle velocity 
				currentParticle.vy = 0 - currentParticle.vy;
				//put the particle on the ground
				currentParticle.y = earthDisplacement + .00001 + currentParticle.radius;				
				}
			}

			//we see what we have to do for collisions with the box.
			//does not take care of collisions from above the box correctly yet
			if(boxOn){

				double zpos = currentParticle.z;
				double xpos = currentParticle.x;
				double r = currentParticle.radius;

				double check = zpos + r - boxSize;
				if((check >= 0 && check <= 2*r) && xpos <= boxSize && xpos >= -boxSize){
					currentParticle.vz = -currentParticle.vz;
				}

				check = zpos - r + boxSize;
				if((check <= 0 && check >= -2*r) && xpos <= boxSize && xpos >= -boxSize){
					currentParticle.vz = -currentParticle.vz;
				}

				check = xpos - r + boxSize;
				if((check <= 0 && check >= -2*r) && zpos <= boxSize && zpos >= -boxSize){
					currentParticle.vx = -currentParticle.vx;
				}

				check = xpos + r - boxSize;
				if((check >= 0 && check <= 2*r) && zpos <= boxSize && zpos >= -boxSize){
					currentParticle.vx = -currentParticle.vx;
				}
			}

			if(bBall){
				double xpos = currentParticle.x;
				double ypos = currentParticle.y;
				double zpos = currentParticle.z;
				double r = currentParticle.radius;

				//first we see if it is in the plane of the backboard
				if (abs(xpos - xDistance) < r){
					//it is in the plane, so now we see if it is in the part of the plane that the board is
					if(abs(zpos) < boardSize){
						if(ypos - boardHeight > 0 && ypos - boardHeight - boardSize <0){
							//the .8 is for damping
							currentParticle.vx = -currentParticle.vx *.8;
						}
					}
					
				}

				//now we have to detect a collision with the torus
				//first we look at its shadow on the xz plane
				//if this shadow intersects the shadow of the torus, there might be a collision

				//we want a frame such that the middle of the torus is at x=0 z=0. The Z part is already done, this is the transform for the X part
				xpos = xpos - xDistance + hoopRadius;

				double shadowRadius = sqrt( xpos*xpos + zpos*zpos);

				if (shadowRadius < (hoopRadius + hoopThickness + r) && shadowRadius > (hoopRadius - hoopThickness - r)){
					//so it is in the shadow of the torus. Next, we find the center of the ring that is closest to the ball
					//that is, if we made a plane through the center of the ball and the middle of the torus, the point we are looking for is the middle of the ring we just cut


					double centerX = xpos/shadowRadius;
					double centerZ = zpos/shadowRadius;
					double centerY = boardHeight;

					double xsquare = (xpos - centerX)*(xpos - centerX);
					double zsquare = (zpos - centerZ)*(zpos - centerZ);
					double ysquare = (ypos - centerY)*(ypos - centerY);

					double distance = sqrt(xsquare + ysquare + zsquare);

					if(distance < hoopThickness + r){
						//currentParticle.vz = -currentParticle.vx;

						//first, we get a vector of the direction pointing out from the ring to the ball
						double vecx = -centerX + xpos;
						double vecy = -centerY + ypos;
						double vecz = -centerZ + zpos;
						double veclength = sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
						
						//now we unit length the vector
						vecx = vecx/veclength;
						vecy = vecy/veclength;
						vecz = vecz/veclength;

						if(vecx < REALLY_SMALL) vecx = 0;
						if(vecy < REALLY_SMALL) vecy = 0;
						if(vecz < REALLY_SMALL) vecz = 0;

						//get the components of the speed
						double tempvx = currentParticle.vx;
						double tempvy = currentParticle.vy;
						double tempvz = currentParticle.vz;

						//double totalv = sqrt(tempvx*tempvx + tempvy*tempvy + tempvz*tempvz);

						//unit length them
						//tempvx = tempvx/totalv;
						//tempvy = tempvy/totalv;
						//tempvz = tempvz/totalv;

						//now we dot product the velocity components and the unit vector
						double dotx = tempvx*vecx;
						double doty = tempvy*vecy;
						double dotz = tempvz*vecz;

						double dotProd = dotx + doty + dotz;

					
						//damping is how much bounce we want off the rim. 2 is a perfect reflection, over 2 adds energy, unfer 2 takes it away. 1 or less will act weird
						double damping = 2.0;
						currentParticle.vx -= damping*dotProd*vecx;
						currentParticle.vy -= damping*dotProd*vecy;
						currentParticle.vz -= damping*dotProd*vecz;


						//for now, going to freeze the balls
						//currentParticle.vz = 0;
						//currentParticle.vx = 0;
						//currentParticle.vy = 0;
						//currentParticle.mass = 0;
					}
				}
			}
		
		}
		particleWorld[i] = currentParticle;
		
	}


	
}


//if collisions are turned on, we both detect and resolve collisions with this function
Particle collision(Particle currentParticle, Particle compParticle, double distanceBetween, int i, int j){
	if((distanceBetween < currentParticle.radius + compParticle.radius) && (i < j) ){
		//if the distance between is less then the combined radius' of the 2 particles, a collison happens
		// i<j is so this is only calculated once

		//now that we know the collision is going to happen, we need to know the velocities of
		//the 2 particles before the collision

		//First we get all the info
		double v1x = currentParticle.vx;
		double v1y = currentParticle.vy;
		double v1z = currentParticle.vz;
		double v2x = compParticle.vx;
		double v2y = compParticle.vy;
		double v2z = compParticle.vz;
		double m1 = currentParticle.mass;
		double m2 = compParticle.mass;


		//now we compute the vector connecting the 2 centers of the spheres. this will determine how the collision happens
		//u is the vecotr connecting current particle to comp particle
		// -u describes comp to current
		double u[3];
		u[0] = compParticle.x - currentParticle.x;
		u[1] = compParticle.y - currentParticle.y;
		u[2] = compParticle.z - currentParticle.z;

		//time to find the unit vector
		double u_length = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
		double unit_u[3];
		unit_u[0] = u[0]/u_length;
		unit_u[1] = u[1]/u_length;
		unit_u[2] = u[2]/u_length;

		//now we dot u with current particles velocity to get the component we will be "swapping" in the collision
		double swap_v_size1 = unit_u[0]*v1x + unit_u[1]*v1y + unit_u[2]*v1z;
		//we have the magnitude of the swap, and the direction will be in direction u, so lets put 1 and 1 together
		double swap_v1[3];
		//this if statement is to stop roundoff errors that occur in glancing collisions
		if(swap_v_size1 > REALLY_SMALL){
			swap_v1[0] = unit_u[0]*swap_v_size1;
			swap_v1[1] = unit_u[1]*swap_v_size1;
			swap_v1[2] = unit_u[2]*swap_v_size1;
		}else{
			swap_v1[0] = 0;
			swap_v1[1] = 0;
			swap_v1[2] = 0;
		}
		//we subtract this swap velocity from the total velocity to get the "remaining" velocity
		double current_v[3];
		current_v[0] = v1x - swap_v1[0];
		current_v[1] = v1y - swap_v1[1];
		current_v[2] = v1z - swap_v1[2];


		//do it all again for comp particle mostly the same, but u is pointing in the opposite direction
		u[0] = -compParticle.x + currentParticle.x;
		u[1] = -compParticle.y + currentParticle.y;
		u[2] = -compParticle.z + currentParticle.z;

		//unit length
		u_length = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
		unit_u[0] = u[0]/u_length;
		unit_u[1] = u[1]/u_length;
		unit_u[2] = u[2]/u_length;

		//dot product
		double swap_v_size2 = (unit_u[0]*v2x + unit_u[1]*v2y + unit_u[2]*v2z);
		double swap_v2[3];
		if(swap_v_size2 > REALLY_SMALL){
			swap_v2[0] = unit_u[0]*swap_v_size2;
			swap_v2[1] = unit_u[1]*swap_v_size2;
			swap_v2[2] = unit_u[2]*swap_v_size2;
		}else{
			swap_v2[0] = 0;
			swap_v2[1] = 0;
			swap_v2[2] = 0;
		}
		//subtract the swap component off
		double comp_v[3];
		comp_v[0] = v2x - swap_v2[0];
		comp_v[1] = v2y - swap_v2[1];
		comp_v[2] = v2z - swap_v2[2];




		//we now do the swap, and turn the swapped momentums into velocities, then add it all together
		current_v[0] = current_v[0] + swap_v2[0]*m2/m1;
		current_v[1] = current_v[1] + swap_v2[1]*m2/m1;
		current_v[2] = current_v[2] + swap_v2[2]*m2/m1;

		comp_v[0] = comp_v[0] + swap_v1[0]*m1/m2;
		comp_v[1] = comp_v[1] + swap_v1[1]*m1/m2;
		comp_v[2] = comp_v[2] + swap_v1[2]*m1/m2;

		//now we change the values accordingly
		currentParticle.vx = current_v[0];
		currentParticle.vy = current_v[1];
		currentParticle.vz = current_v[2];
		compParticle.vx = comp_v[0];
		compParticle.vy = comp_v[1];
		compParticle.vz = comp_v[2];

		//save the info to particle world for particle j (comp particle)
		particleWorld[j] = compParticle;
	} 
	return currentParticle;
}

//deals with key hits
void keyhit(unsigned char key, int x, int y)
{
	Color white;
	white.red = 1;
	white.green = 1;
	white.blue = 1;

	Color orange;
	orange.red = 1;
	orange.green = .55;
	orange.blue = 0;

	switch(key) {
	case 'b':
		boxOn = !boxOn;
		break;
		//wooo, basketball time
	case 'B':
		bBall = !bBall;
		break;
	case 'k':
		collisionsOn = !collisionsOn;
		break;
	case ' ':
		animationOn = !animationOn;
		break;
		//forward
	case 'w':{
		volatile float mouse1Rad = mouse1*PI/180;
		volatile float mouse2Rad = mouse2*PI/180;

		volatile float xRotation1 = sin(mouse1Rad);
		volatile float xRotation2 = cos(mouse1Rad);
		
		volatile float yRotation1 = sin(mouse2Rad);
		volatile float yRotation2 = cos(mouse2Rad);


		worldX += .25*xRotation2;
		worldZ += .25*xRotation1;
		break;
			 }
		//backward
	case 's':
		{
		volatile float mouse1Rad = mouse1*PI/180;
		volatile float mouse2Rad = mouse2*PI/180;

		volatile float xRotation1 = sin(mouse1Rad);
		volatile float xRotation2 = cos(mouse1Rad);
		
		volatile float yRotation1 = sin(mouse2Rad);
		volatile float yRotation2 = cos(mouse2Rad);


		worldX += -.25*xRotation2;
		worldZ += -.25*xRotation1;
		break;
			 }
		//left
	case 'a':
		{
		volatile float mouse1Rad = mouse1*PI/180;
		volatile float mouse2Rad = mouse2*PI/180;

		volatile float xRotation1 = sin(mouse1Rad);
		volatile float xRotation2 = cos(mouse1Rad);
		
		volatile float yRotation1 = sin(mouse2Rad);
		volatile float yRotation2 = cos(mouse2Rad);


		worldZ -= .25*xRotation2;
		worldX += .25*xRotation1;
		break;
			 }
		//right
	case 'd':
		{
		volatile float mouse1Rad = mouse1*PI/180;
		volatile float mouse2Rad = mouse2*PI/180;

		volatile float xRotation1 = sin(mouse1Rad);
		volatile float xRotation2 = cos(mouse1Rad);
		
		volatile float yRotation1 = sin(mouse2Rad);
		volatile float yRotation2 = cos(mouse2Rad);


		worldZ += .25*xRotation2;
		worldX -= .25*xRotation1;
		break;
			 }
		//up
	case 'z':
		worldY += .25;
		break;
		//down
	case 'x':
		worldY -= .25;
		break;
	case 'e':
		earthOn = !earthOn;
		break;
	case 'C':
		particleWorld.clear();
		break;
	case 'p':
		newBall = !newBall;
		break;
	case '>':
		speed ++;
		break;
	case '<':
		speed--;
		break;
	case 'P':{
		//shoot out a new ball
		//it will be shot where the camera is facing
		float vx, vy, vz;
		if (speed == 0){
			vx = 0; 
			vy = 0; 
			vz = 0;
		}else{
			float mouse1Rad = mouse1*PI/180;
			float mouse2Rad = mouse2*PI/180;

			float xRotation1 = sin(mouse1Rad);
			float xRotation2 = cos(mouse1Rad);

			float yRotation1 = sin(mouse2Rad);
			float yRotation2 = cos(mouse2Rad);

			vx = xRotation2;
			vy = -yRotation1;
			vz = (xRotation1);

			//speed affects how fast we shoot it
			float normalizer = sqrt(vx*vx + vy*vy + vz*vz)/speed;

			vx = vx/normalizer;
			vy = vy/normalizer;
			vz = vz/normalizer;
		}
		Particle newParticle = particle(10000000000, .000001, .4, centerX, centerY, centerZ, vx, vy, vz, white);
		break;}
	case '1':{
		particleWorld.clear();
		Particle part1 = particle(1,0,.1,0,-0.3,-.3,0,0,.3,white);
		Particle part2 = particle(1,0,.1,0,-0.3,.3,0,0,-.3,orange);
		Particle part3 = particle(1,0,.1,0.3,-0.3,0,0.3,0,0,orange);
		Particle part4 = particle(1,0,.1,-0.3,-0.3,0,0-.3,0,0,orange);

		break;}
	case '2':{
		particleWorld.clear();
		Particle vertDrop1 = particle(1,0,.1,0,0,0,0,0,0,white);
		Particle vertDrop2 = particle(1,0,.1,0,.6,0,0,0,0,orange);
		break;}
	case '3':{
		particleWorld.clear();
		Particle skim1 = particle(1,0,.1,0,0,0,0,0,0,white);
		Particle skim2 = particle(1,0,.1,.199,.6,0,0,0,0,orange);
		break;}
	}
}


//dealing with mouse clicks
void mouseClick(int button, int state, int x, int y){
	mouseX = x;
	mouseY = y;

}

void motion(int x, int y){
	float changeX = (mouseX - x)/4;
	float changeY = (mouseY - y)/4;
	//divide by 4 just to damp the motion a bit

	mouse1 = mouse1 + changeX;
	mouse2 = mouse2 + changeY;

	mouseX = x;
	mouseY = y;


}

// Draw frame
void DrawWorld()
{
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
	
	//turn our degree mouse turns into radians
	float mouse1Rad = mouse1*PI/180;
	float mouse2Rad = mouse2*PI/180;

	float xRotation1 = sin(mouse1Rad);
	float xRotation2 = cos(mouse1Rad);

	float yRotation1 = sin(mouse2Rad);
	float yRotation2 = cos(mouse2Rad);

	x = -5 + worldX;
	y = 0 + worldY;
	z = 0 + worldZ;
	
	float viewDistance = 5;

	//position of the refrence point
	centerX = viewDistance*(xRotation2) + x;
	centerY = -yRotation1*viewDistance + y;
	centerZ = viewDistance*(xRotation1) + z;
	//position of the eye
	float eye[] = {x,y,z};
	
	//direction of the up vector
	float up[] = {0,1,0};

	//set up the camera where i want with this
	gluLookAt(eye[0],eye[1],eye[2],
			centerX,centerY,centerZ,
			up[0],up[1], up[2]);

	// rotate around the y-axis according to the mouse movement
	//glRotatef(mouse1,0,1,0);
	//now about the x axis
	//glRotatef(mouse2, 1, 0, 0);
	//glTranslatef(worldX,0,worldZ);	

	// set the background colour to be used, then apply it
	glClearColor(0,0,0,0);
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//glTranslatef(worldX,0,worldZ);

	//drawing the particles here
	for (int i = 0 ; i < particleWorld.size() ; i++){
		Particle currentParticle = particleWorld[i];
		Color partColor = currentParticle.color;

		glPushMatrix();
		glColor3f(partColor.red,partColor.green,partColor.blue);
		glTranslatef(currentParticle.x,currentParticle.y,currentParticle.z);
		glutSolidSphere(currentParticle.radius,sphereQuality,sphereQuality);
		glPopMatrix();
		
	}

	//draw the ghost ball if we are in new ball mode
	if (newBall){
		glPushMatrix();
		glColor3f(0.5,0.5,0.5);
		glTranslatef(centerX,centerY,centerZ);
		glutSolidSphere(0.1,sphereQuality,sphereQuality);
		glPopMatrix();
	}
	glPushMatrix();
	//glLoadIdentity();

	// draw the earth if its on
	if(earthOn){
		glBegin(GL_QUADS);
		glColor3f(1,0,0);
		glNormal3f(0,1.0,0);
		glVertex3f( -10.0f, earthDisplacement, -10.0f );

		glColor3f(0,1,0);
		glNormal3f(0,1.0,0);
		glVertex3f( 10.0f, earthDisplacement, -10.0f );

		glColor3f(0,0,1);
		glNormal3f(0,1.0,0);
		glVertex3f( 10.0f, earthDisplacement, 10.0f );

		glColor3f(1,1,1);
		glNormal3f(0,1.0,0);
		glVertex3f( -10.0f, earthDisplacement, 10.0f );
		glEnd();
	}

	if (boxOn){
		
		

		glBegin(GL_QUADS);
		glColor3f(1,0,0);

		glVertex3f( boxSize, earthDisplacement, boxSize );

		glColor3f(0,1,0);
		glVertex3f( -boxSize, earthDisplacement, boxSize );

		glColor3f(0,0,1);
		glVertex3f( -boxSize, boxHeight, boxSize );

		glColor3f(1,1,1);
		glVertex3f( boxSize, boxHeight, boxSize );
		glEnd();

		glBegin(GL_QUADS);
		glColor3f(1,0,0);
		glVertex3f( -boxSize, earthDisplacement, boxSize );

		glColor3f(0,1,0);
		glVertex3f( -boxSize, earthDisplacement, -boxSize );

		glColor3f(0,0,1);
		glVertex3f( -boxSize, boxHeight, -boxSize );

		glColor3f(1,1,1);
		glVertex3f( -boxSize, boxHeight, boxSize );
		glEnd();

		glBegin(GL_QUADS);
		glColor3f(1,0,0);
		glVertex3f( boxSize, earthDisplacement, -boxSize );

		glColor3f(0,1,0);
		glVertex3f( boxSize, earthDisplacement, boxSize );

		glColor3f(0,0,1);
		glVertex3f( boxSize, boxHeight, boxSize );

		glColor3f(1,1,1);
		glVertex3f( boxSize, boxHeight, -boxSize );
		glEnd();

		glBegin(GL_QUADS);
		glColor3f(1,0,0);
		glVertex3f( boxSize, earthDisplacement, -boxSize );

		glColor3f(0,1,0);
		glVertex3f( -boxSize, earthDisplacement, -boxSize );

		glColor3f(0,0,1);
		glVertex3f( -boxSize, boxHeight, -boxSize );

		glColor3f(1,1,1);
		glVertex3f( boxSize, boxHeight, -boxSize );
		glEnd();

	}

	if(bBall){
		//double xDistance = 10;
		//double boardSize = 3;
		//double hoopRadius = 1;
		//double boardHeight =3;

		//the backboard
		glBegin(GL_QUADS);
		glColor3f(.93, .863, .508);
		//glColor3f(1, 1, 1);
		glVertex3f( xDistance, boardHeight, -boardSize/2 );

		glVertex3f( xDistance, boardHeight, boardSize/2 );

		glVertex3f( xDistance, boardHeight + boardSize, boardSize/2 );

		glVertex3f( xDistance, boardHeight + boardSize, -boardSize/2 );
		glEnd();

		//the hoop
		glPushMatrix();
		glColor3f(0.8,0.8,0.8);
		
		glTranslatef(xDistance - hoopRadius,boardHeight,0);
		glRotatef(90, 1.0, 0, 0);
		double thickness = .2;
		glutSolidTorus(hoopThickness, hoopRadius, 20, 20);
		glPopMatrix();

	}

	//glColor3f(0.8, 0.8,0.8);      // choose a light-grey color
	//glTranslatef(0,0,0.0);        // decided not to move the torus
	//glutSolidTorus(0.1,0.5,8,20); // draw a solid torus



	glPopMatrix();

/*	GLfloat light_position[] = { 0.0, 1.0, 0.0, 0.0 };
	GLfloat ambient_light[] = {.2, .2, .2, 1.0};
	GLfloat diffuse_light[] = {.8,.8,.8,1.0};
	glShadeModel (GL_SMOOTH);
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient_light);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse_light);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_DEPTH_TEST);
	*/
	glutSwapBuffers();          // force the draw-buffer to be displayed
}


void myReshape(int w, int h)
{
	float Fovy = 35.0;    // field of view  in the y-direction (in degrees)
	float Znear = 0.2;    // near clipping plane (distance from camera)
	float Zfar = 20.0;   // far clipping plane (distance from camera)
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(Fovy, (GLfloat) w/h, Znear, Zfar);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

//main loop
//doing a few things im not sure about, but it seems to work
int WINAPI WinMain (HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow)
{


	//this should NEED to be called, but when we do, it breaks it. somehow works without an initializer being called
	//glutInit( &argc, argv );
	glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
	glutInitWindowSize( 640, 480 );
	mainWindow = glutCreateWindow( "CJ's Crazy Awesome Physics Simulator" );
	glutReshapeFunc(myReshape);
	glutDisplayFunc( DrawWorld );
	glutIdleFunc(Idle);
	glutKeyboardFunc(keyhit);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_COLOR_MATERIAL);

	//get some basic lighting going
	GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat mat_shininess[] = { 50.0 };
	GLfloat light_position[] = { 100.0, 100.0, 10.0, 0.0 };
	glClearColor (0.0, 0.0, 0.0, 0.0);
	glShadeModel (GL_SMOOTH);

	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);




	//added these for mouse functions
	glutMouseFunc(mouseClick);
	glutMotionFunc(motion);

	
/*	trying to set up glui stuff, not going too well

	//  pointer to a GLUI window
	GLUI * glui_window;
	//  Create GLUI window
	glui_window = GLUI_Master.create_glui ("Options", GLUI_DOUBLE , 0, 0);

	glui_window->set_main_gfx_window (mainWindow);
	GLUI_Master.set_glutIdleFunc (Idle);

*/	
	// pass control over to GLUT
	glutMainLoop();

	return 0;       // never reached
}



