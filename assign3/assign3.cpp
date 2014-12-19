/*
CSCI 480
Assignment 3 Raytracer

Name: <Your name here>
*/

#include <pic.h>
#include <windows.h>
#include <stdlib.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <vector>
#include <stdio.h>
#include <string>
#include <iostream>

#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

char *filename=0;
//some methods declared
bool collision(double origin[3],double point[3],std::string &OBJ,int &index,double *iPoint, std::string currObj, int currIndex);
std::vector<double> phongTri(double p[3],int index);

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//resolution different than width and height for debugging purposes
#define resW 640
#define resH 480

//the field of view of the camera
#define fov 60.0
#define PI 3.141592
  //Tan takes radians so we need to convert fov to radfov
  double radFov = fov*(double)2*PI/(double)360;

//the aspect ratio
float alpha = (float)WIDTH/HEIGHT;

unsigned char buffer[HEIGHT][WIDTH][3];

//reflection constant
#define rconst 2

struct Vertex
{
  double position[3];//x,y,z
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
  bool hit; 
};


//2D Matrix of Image Plane Pixels
Vertex pix[resW][resH];

typedef struct _Triangle
{
  struct Vertex v[3];
} Triangle;

typedef struct _Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
} Sphere;

typedef struct _Light
{
  double position[3];
  double color[3];
} Light;

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles=0;
int num_spheres=0;
int num_lights=0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

//origin-camera position
double cam[3] = {0,0,0};//x y z


//Normalize a vector
void normalize(double *n){
	 
	 double size = sqrt( (n[0] * n[0]) + (n[1] * n[1]) + (n[2] * n[2]) );

		  //normalized
		  n[0] /= size;
		  n[1] /= size;
		  n[2] /= size;


 }

//Dot product of two vectors
double dot( double x[3] ,double y[3]){
	
	  double answer =   (x[0]*y[0]) + (x[1]*y[1]) + (x[2]*y[2]);
	   
	  return answer;

}

//Cross product of two vectors
void cross(double  *n, double a[3], double b[3]){
	
	n[0] = ((a[1]*b[2]) - (a[2]*b[1]));
	n[1] = (a[2]*b[0]) - (a[0]*b[2]);
	n[2] =(a[0]*b[1]) - (a[1]*b[0]);


}

//finds the direction vector given two points 
void direction(double start[3], double end[3], double* dir){
	
	dir[0] = end[0] - start[0];
	dir[1] = end[1] - start[1];
	dir[2] = end[2] - start[2];

}

 bool check_same_clock_dir(double pt1[3], double pt2[3],  double pt3[3], double norm[3])
{  
   float testi, testj, testk;
   float dotprod;
   // normal of trinagle
   testi = (((pt2[1] - pt1[1])*(pt3[2] - pt1[2])) - ((pt3[1] - pt1[1])*(pt2[2] - pt1[2])));
   testj = (((pt2[2] - pt1[2])*(pt3[0] - pt1[0])) - ((pt3[2] - pt1[2])*(pt2[0] - pt1[0])));
   testk = (((pt2[0] - pt1[0])*(pt3[1] - pt1[1])) - ((pt3[0] - pt1[0])*(pt2[1] - pt1[1])));

   // Dot product with triangle normal
   dotprod = testi*norm[0] + testj*norm[1] + testk*norm[2];
   
   //answer
   if(dotprod < 0) return false;
   else return true;
}

//ray-triangle intersection
bool iTriangle(double o[3], double p[3], Triangle tri, double *iPoint){ //given the line endpoints, and triangle we can find if a ray intercets with a triangle

	
	//1. get triangle normal by crossing two edges, and normalize it
	double first[3];
	double second[3];
	direction(tri.v[0].position, tri.v[1].position,first);
	direction(tri.v[0].position, tri.v[2].position,second);

	double triNorm[3];
    cross(triNorm,first,second);
	normalize(triNorm);

	//2.get direction of ray and normalize it
	double dir[3];
	direction(o,p,dir);
	normalize(dir);
		

	//3.check if ray intersects plane

	double t = dot(triNorm,dir);
	
	//3a. plane-ray intersection
	if(t == 0) 			
		return false;
	
	//3b. find intersection point

	double oma[3];
	direction(tri.v[0].position,o,oma);
	t =  -1*dot(triNorm,oma);
	t/= dot(triNorm,dir);
		
	if(t<0)
		return false; //no intersection

	
	iPoint[0] = o[0]+(t*dir[0]);
	iPoint[1]=  o[1]+(t*dir[1]);
    iPoint[2]=  o[2]+(t*dir[2]);

	//3c.check if intersection point lies inside triangle using clockwise rule - found online
		if(check_same_clock_dir(tri.v[0].position, tri.v[1].position, iPoint, triNorm) == true)
      {
		  if(check_same_clock_dir(tri.v[1].position, tri.v[2].position, iPoint, triNorm) == true)
         {
			 if(check_same_clock_dir(tri.v[2].position, tri.v[0].position, iPoint, triNorm) == true)
            {   
				
               return true; //ray hits inside triangle
            }
         }
      }

	return false; //ray intersection lies outside triangle
}

//ray-sphere intersection
bool iSphere(double o[3], double p[3], Sphere s, double *iPoint){ //given 2 endpoints of line, and a sphere to check intersetions
	 double t1,t2;

	//1.calculate the direction of the viewRay
	 double currViewRay[3];
	 direction(p,o,currViewRay);
	 normalize(currViewRay);

    //1a. caluclate direction of vector pointing to center
	double v[3];
	direction(s.position,o,v);

	
	//2. calculate a, b, c to plug into "b^2-4ac"
	double A,B,C; 
	A = dot( currViewRay, currViewRay);
	B =  2 * dot(currViewRay,v);
	C  = dot(v,v) - (s.radius * s.radius);

	if((B*B - 4*A*C) < 0)
		return false; //no ray-sphere intersection
	else{
		//calculating the intersection point
		//if discriminat is zero there is one intersection
		//else there are two
		//found by t = (-b +/- sqrt(b2 - 4 a c)) / (2 a)

		if((B*B - 4*A*C) == 0){
			t1 = (-B)/(2*A);
		    iPoint[0]= o[0]+(t1*currViewRay[0]);
			iPoint[1]= o[1]+(t1*currViewRay[1]);
			iPoint[2]= o[2]+(t1*currViewRay[2]);
			return true;
          
		}

		else{
		    //note: t1 is closer than t2
           t1=(.5)*((-B)- sqrt(B*B - 4*A*C))/(A);
		   t2=.5 * ((-B)+ sqrt(B*B - 4*A*C))/(A);

		   //t1<0: intersects behind camera, so t2 is the answer
		   if(t1<0){
            iPoint[0]= o[0]+(t2*currViewRay[0]);
			iPoint[1]= o[1]+(t2*currViewRay[1]);
			iPoint[2]= o[2]+(t2*currViewRay[2]);
			return true; //ray-sphere intersection found!

		   }
		   //t1: is infron of camera, therefore the answer
		   else{
            iPoint[0]= o[0]+(t1*currViewRay[0]);
			iPoint[1]= o[1]+(t1*currViewRay[1]);
			iPoint[2]= o[2]+(t1*currViewRay[2]);
            return true; //ray-sphere intersection found!
		   }


        

		}

	}
}

//check for any shadowRay intersection to a light source
bool shadowCheck(double iPoint[3], double j, std::string o,Light sun){
     //check intersection of shadworay- with every other object for every light source
	double check[3];
    double light[3];
	double behind[3];
	double s1;
	double s2;
	bool shade = false;

	  
	     

			       //check shadowRay-triangle intersection between ipoint and any lightsource
		   for(int k = 0; k<num_triangles; k++){	  
				 if(iTriangle(iPoint,sun.position,triangles[k],check)){
                      if(o != "t" || j!=k){
					  direction(sun.position,iPoint,light); 
					  direction(sun.position,check,behind);  
					  s1 = sqrt( (light[0] * light[0]) + (light[1] * light[1]) + (light[2] * light[2]) );
					  s2 =  sqrt( (behind[0] * behind[0]) + (behind[1] * behind[1]) + (behind[2] * behind[2]) );
					  if(s1>s2)
					 return true;
					  }
			    }	
		   }
		
		  //checks shadowRay-sphere intersection
		 if(!shade){//havent found shade yet
		  for(int k = 0; k<num_spheres; k++){
		        if(iSphere(iPoint,sun.position,spheres[k],check)){
				   if(o != "s" || j!=k){
					 direction(sun.position,iPoint,light); 
					  direction(sun.position,check,behind);  
					  s1 = sqrt( (light[0] * light[0]) + (light[1] * light[1]) + (light[2] * light[2]) );
					  s2 =  sqrt( (behind[0] * behind[0]) + (behind[1] * behind[1]) + (behind[2] * behind[2]) );
					  if(s1>s2){
					   return true;
					  }
				   }
				}
		  }
		 }
		
		   //if no shade is found, the pixel has atleast one light source hitting it, so it will be not under a shadow
		 

	return false;
}

//phong model to illuminate a sphere
std::vector<double> phongSphere(double p[3], int index){
	std::vector<double> Color(3);
	double r[3];
	double I[3]={0,0,0};
	double l[3]={0,0,0};
	 //1.find surface normal "n"
	double n[3];
	direction(spheres[index].position,p,n);
	normalize(n);

	//2.find viewRay
    double v[3];
	direction(p,cam,v);
	normalize(v);

    //3.go thru all existing light sources
	for(int f=0;f<num_lights;f++){



    if(!shadowCheck(p,index,"s",lights[f])){

    //4.find lightVector "l"
	direction(p,lights[f].position,l);
	normalize(l);
	
	//3.find ray vector "r"
	r[0]= 2*dot(l,n)*n[0]-l[0];
	r[1]= 2*dot(l,n)*n[1]-l[1];
	r[2]= 2*dot(l,n)*n[2]-l[2];
	normalize(r);
	//compute phong illumination equation
	double D[3];//diffused
	double S[3];//specular
	
	double lambert = dot(n,l);
	if(lambert<0)
		lambert =0;

	//diffused lighting: difused_color + lambert coef
	D[0]=  spheres[index].color_diffuse[0] * lambert;
	D[1]=  spheres[index].color_diffuse[1] * lambert;
	D[2]=  spheres[index].color_diffuse[2] * lambert;

	double rv = dot(r,v);
	if(rv<0)
     rv=0;

	//specular lighting: specular_color * dot(raycastVector, viewRay)
	S[0]= spheres[index].color_specular[0] * pow(rv,spheres[index].shininess);
	S[1]= spheres[index].color_specular[1] * pow(rv,spheres[index].shininess);
	S[2]= spheres[index].color_specular[2] * pow(rv,spheres[index].shininess);

	//for each lightsource, we compute I for an added effect
	
	I[0] += lights[f].color[0] * (D[0] + S[0]);
    I[1] += lights[f].color[1] * (D[1] + S[1]);
	I[2] += lights[f].color[2] * (D[2] + S[2]);
	}
	}
	//add up calculations to determine final color
	 Color[0]= (ambient_light[0] + I[0]);
	 Color[1]= (ambient_light[1] + I[1]);
	 Color[2]= (ambient_light[2] + I[2]);



	//else the object is assumed to be non-reflective
	

	return Color;



}

double area(double a[3], double b[3], double c[3]){
	double side1[3];
	double side2[3];
	double side3[3];
	//sizes of each side
	double sizeA,sizeB,sizeC; 
	
	direction(a,b,side1);
	direction(a,c,side2);
	direction(c,b,side3);

	sizeA = sqrt( (side1[0] * side1[0]) + (side1[1] * side1[1]) + (side1[2] * side1[2]) );
	sizeB = sqrt( (side2[0] * side2[0]) + (side2[1] * side2[1]) + (side2[2] * side2[2]) );
	sizeC = sqrt( (side3[0] * side3[0]) + (side3[1] * side3[1]) + (side3[2] * side3[2]) );

	double s = (.5) * (sizeA + sizeB + sizeC);


	return sqrt(s * (s-sizeA) * (s-sizeB) * (s-sizeC));

}

std::vector<double> phongTri(double p[3],int index){
	std::vector<double> Color(3);
	double I[3]={0,0,0};
	double l[3]={0,0,0};
	//compute phong illumination equation vars
	double D[3];//diffused
	 double d1[3];
     double d2[3];
	 double d3[3];
	double S[3];//specular
	 double s1[3];
	 double s2[3];
	 double s3[3];
	 double r[3];
	 double lambert, rv;

	//1.find viewRay
    double v[3];
	direction(p,cam,v);
	normalize(v);

	//2. find areas
	double A1 =	area(triangles[index].v[1].position, triangles[index].v[2].position,p); //points: 1,2,p;
	double A2 = area(triangles[index].v[0].position, triangles[index].v[2].position,p);//points: 0,2,p
	double A3 =  area(triangles[index].v[0].position, triangles[index].v[1].position,p);//points: 0,1,p
	double totsA = 	area(triangles[index].v[0].position, triangles[index].v[1].position,triangles[index].v[2].position); //points: 0,1,2

	//3.find surface normal "n" of a triangle(linear intropolation)
	double n1[3];
	double n2[3];
	double n3[3];
	double n[3];
	//linear intropolate 3 normals
	n[0] = ((A1/totsA)*triangles[index].v[0].normal[0]) + ((A2/totsA)*triangles[index].v[1].normal[0]) + ((A3/totsA)*triangles[index].v[2].normal[0]);
	n[1] = ((A1/totsA)*triangles[index].v[0].normal[1]) + ((A2/totsA)*triangles[index].v[1].normal[1]) + ((A3/totsA)*triangles[index].v[2].normal[1]);
	n[2] = ((A1/totsA)*triangles[index].v[0].normal[2]) + ((A2/totsA)*triangles[index].v[1].normal[2]) + ((A3/totsA)*triangles[index].v[2].normal[2]);



	normalize(n);
	



	//for each light source..
    //4.find lightVector "l"

    for(int f=0; f<num_lights;f++){

		//check if it hits shadow
		if(!shadowCheck(p,index,"t",lights[f])){


		
	direction(p,lights[f].position,l);
	normalize(l);

	//3.find ray vector "r"
	 r[0]= 2*dot(l,n)*n[0]-l[0];
	 r[1]= 2*dot(l,n)*n[1]-l[1];
     r[2]= 2*dot(l,n)*n[2]-l[2];
	
	
	lambert = dot(n,l);
	if(lambert<0)
		lambert =0;

	//diffused lighting: difused_color + lambert coef (for each tri vertex)
	d1[0]= triangles[index].v[0].color_diffuse[0] * lambert;
	d1[1]= triangles[index].v[0].color_diffuse[1] * lambert;
	d1[2]= triangles[index].v[0].color_diffuse[2] * lambert;

	d2[0]= triangles[index].v[1].color_diffuse[0] * lambert;
	d2[1]= triangles[index].v[1].color_diffuse[1] * lambert;
	d2[2]= triangles[index].v[1].color_diffuse[2] * lambert;


	d3[0]= triangles[index].v[2].color_diffuse[0] * lambert;
	d3[1]= triangles[index].v[2].color_diffuse[1] * lambert;
	d3[2]= triangles[index].v[2].color_diffuse[2] * lambert;


	 rv = dot(r,v);
	if(rv<0)
     rv=0;

	//specular lighting: specular_color * dot(raycastVector, viewRay)
	s1[0]= triangles[index].v[0].color_specular[0] * pow(rv,triangles[index].v[0].shininess);
	s1[1]= triangles[index].v[0].color_specular[1] * pow(rv,triangles[index].v[0].shininess);
	s1[2]= triangles[index].v[0].color_specular[2] * pow(rv,triangles[index].v[0].shininess);

	s2[0]= triangles[index].v[1].color_specular[0] * pow(rv,triangles[index].v[1].shininess);
	s2[1]= triangles[index].v[1].color_specular[1] * pow(rv,triangles[index].v[1].shininess);
	s2[2]= triangles[index].v[1].color_specular[2] * pow(rv,triangles[index].v[1].shininess);

    s3[0]= triangles[index].v[2].color_specular[0] * pow(rv,triangles[index].v[2].shininess);
	s3[1]= triangles[index].v[2].color_specular[1] * pow(rv,triangles[index].v[2].shininess);
	s3[2]= triangles[index].v[2].color_specular[2] * pow(rv,triangles[index].v[2].shininess);

	//interpolate the values
	D[0] = ((A1/totsA)*d1[0]) + ((A2/totsA)*d2[0]) + ((A3/totsA)*d3[0]);
	D[1] = ((A1/totsA)*d1[1]) + ((A2/totsA)*d2[1]) + ((A3/totsA)*d3[1]);
	D[2] = ((A1/totsA)*d1[2]) + ((A2/totsA)*d2[2]) + ((A3/totsA)*d3[2]);

	S[0] = ((A1/totsA)*s1[0]) + ((A2/totsA)*s2[0]) + ((A3/totsA)*s3[0]);
	S[1] = ((A1/totsA)*s1[1]) + ((A2/totsA)*s2[1]) + ((A3/totsA)*s3[1]);
	S[2] = ((A1/totsA)*s1[2]) + ((A2/totsA)*s2[2]) + ((A3/totsA)*s3[2]);


	I[0] += lights[f].color[0] * (D[0] + S[0]);
    I[1] += lights[f].color[1] * (D[1] + S[1]);
	I[2] += lights[f].color[2] * (D[2] + S[2]);
		}
	}
	//add up calculations to determine final color
	
	Color[0]= (ambient_light[0] + I[0]);
	Color[1]= (ambient_light[1] + I[1]);
	Color[2]= (ambient_light[2] + I[2]);

	/* MY ATTEMPT AT REFLECTION
	if(triangles[index].v[0].shininess >50) //we assume that object is reflective
	{//checks if there is anything to reflect
		double temp[3]; //second point
		temp[0] = p[0] + (1 * r[0]);
		temp[1] = p[1] + (1 * r[1]);
		temp[2] = p[2] + (1 * r[2]);
		double newP[3];//a new intersection point
	    std::string obj;
		int inx;
		std::vector<double> reflected(3);

		if(collision(p,temp,obj,inx,newP,"t",index)){
			
			if(obj=="s"){
				reflected = phongSphere(newP,inx);
				 Color[0]= ((1-triangles[index].v[0].color_specular[0])*Color[0]) + ((triangles[index].v[0].color_specular[0])*reflected[0]);
	             Color[1]= ((1-triangles[index].v[0].color_specular[1])*Color[1]) + ((triangles[index].v[0].color_specular[1])*reflected[1]);
	             Color[2]= ((1-triangles[index].v[0].color_specular[2])*Color[2]) + ((triangles[index].v[0].color_specular[2])*reflected[2]);

			}
			if(obj=="t"){
				reflected = phongTri(newP,inx);
				 Color[0]= ((1-triangles[index].v[0].color_specular[0])*Color[0]) + ((triangles[index].v[0].color_specular[0])*reflected[0]);
	             Color[1]= ((1-triangles[index].v[0].color_specular[1])*Color[1]) + ((triangles[index].v[0].color_specular[1])*reflected[1]);
	             Color[2]= ((1-triangles[index].v[0].color_specular[2])*Color[2]) + ((triangles[index].v[0].color_specular[2])*reflected[2]);
			
			

			}
		}
	 }*/


	return Color;


}
bool collision(double origin[3],double point[3],std::string &OBJ,int &index,double *iPoint, std::string currObj, int currIndex){
	 bool hit = false;
	 double check[3];
	 double minDist=-1000;
		  //check ray-triangle intersection
		  for(int k = 0; k<num_triangles; k++){	  
			  if(!(currObj == "t" && currIndex==k)){
		       if(iTriangle(origin,point,triangles[k],check)){
				//the ray intersects, and info is updated if intersection is closer to camera	
				if(minDist<check[2]){
					iPoint[0] = check[0];
					iPoint[1] = check[1];
					iPoint[2] = check[2];
					minDist = check[2]; 
					OBJ = "t";
					index=k;
					hit=true;
				}
		      }
			}
		   }

		  //checks ray-sphere intersection
		  for(int k = 0; k<num_spheres; k++){
		    if(iSphere(origin,point,spheres[k],check)){
				//the ray intersects, and info is updated if intersection is closer to camera
				  if(!(currObj == "s" && currIndex==k)){
				if(minDist<spheres[k].position[2]){
					iPoint[0] = check[0];
					iPoint[1] = check[1];
					iPoint[2] = check[2];
		            minDist = spheres[k].position[2];
					OBJ = "s";
					index =k;
					hit=true;
				}	
			  }
		    }
			  
		  }
      
      return hit;
}

void draw_scene()
{
  
 //creates a viewRay for each pixel in that direction 
 //and uses that viewRay to check intercetions with all objects
 //if: it intersects, we will use phong shading for that pixel
 //unless its in shadow
 //else: it will be bgcolor


	std::string OBJ = "object";
    int index = 0;
	double iPoint[3];
	double check [3];

     for(int i = 0; i< resW; i++){
	    glPointSize(2.0);
        glBegin(GL_POINTS);
	   for(int j=0; j<resH; j++){
         //for each pixel: find ray-obj intersection closest to camera
		 
		  pix[i][j].hit=false;

	      pix[i][j].hit = collision(cam,pix[i][j].position,OBJ,index,iPoint,"not aplicable",-1); //checks if the pixel hit anything

		  //ray intersected: pixel will get appropriate color
		    if(pix[i][j].hit){
				//sphere shading
				if(OBJ=="s"){
				        std::vector<double>final(3);
						final = phongSphere(iPoint,index);
                        
						final[0] = min(final[0],1);
						final[1] = min(final[1],1);
						final[2] = min(final[2],1);
						plot_pixel(i,j,final[0]*250,final[1]*250,final[2]*250);



				}

				//triangle shading
				if(OBJ=="t"){	
					  std::vector<double>final(3);
						final = phongTri(iPoint,index);
                        
						final[0] = min(final[0],1);
						final[1] = min(final[1],1);
						final[2] = min(final[2],1);
						plot_pixel(i,j,final[0]*250,final[1]*250,final[2]*250);
					
				}
			  //note: pixel shaded differently according to shape
			}

		    //no intersection: backgorund color
			if(!pix[i][j].hit){
				plot_pixel(i,j,250,250,250);
				
			}

	  }
     glEnd();
     glFlush();
    }//end: pixel loop
   
  printf("Done!\n"); fflush(stdout);
}



void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  buffer[HEIGHT-y-1][x][0]=r;
  buffer[HEIGHT-y-1][x][1]=g;
  buffer[HEIGHT-y-1][x][2]=b;
}

void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
      plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  Pic *in = NULL;

  in = pic_alloc(640, 480, 3, NULL);
  printf("Saving JPEG file: %s\n", filename);

  memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
  if (jpeg_write(filename, in))
    printf("File saved Successfully\n");
  else
    printf("Error in Saving\n");

  pic_free(in);      

}

void parse_check(char *expected,char *found)
{
  if(stricmp(expected,found))
    {
      char error[100];
      printf("Expected '%s ' found '%s '\n",expected,found);
      printf("Parse error, abnormal abortion\n");
      exit(0);
    }

}

void parse_doubles(FILE*file, char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv,"r");
  printf(argv);
  int number_of_objects;
  char type[50];
  int i;
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i",&number_of_objects);

  printf("number of objects: %i\n",number_of_objects);
  char str[200];

  parse_doubles(file,"amb:",ambient_light);

  for(i=0;i < number_of_objects;i++)
    {
      fscanf(file,"%s\n",type);
      printf("%s\n",type);
      if(stricmp(type,"triangle")==0)
	{

	  printf("found triangle\n");
	  int j;

	  for(j=0;j < 3;j++)
	    {
	      parse_doubles(file,"pos:",t.v[j].position);
	      parse_doubles(file,"nor:",t.v[j].normal);
	      parse_doubles(file,"dif:",t.v[j].color_diffuse);
	      parse_doubles(file,"spe:",t.v[j].color_specular);
	      parse_shi(file,&t.v[j].shininess);
	    }

	  if(num_triangles == MAX_TRIANGLES)
	    {
	      printf("too many triangles, you should increase MAX_TRIANGLES!\n");
	      exit(0);
	    }
	  triangles[num_triangles++] = t;
	}
      else if(stricmp(type,"sphere")==0)
	{
	  printf("found sphere\n");

	  parse_doubles(file,"pos:",s.position);
	  parse_rad(file,&s.radius);
	  parse_doubles(file,"dif:",s.color_diffuse);
	  parse_doubles(file,"spe:",s.color_specular);
	  parse_shi(file,&s.shininess);

	  if(num_spheres == MAX_SPHERES)
	    {
	      printf("too many spheres, you should increase MAX_SPHERES!\n");
	      exit(0);
	    }
	  spheres[num_spheres++] = s;
	}
      else if(stricmp(type,"light")==0)
	{
	  printf("found light\n");
	  parse_doubles(file,"pos:",l.position);
	  parse_doubles(file,"col:",l.color);

	  if(num_lights == MAX_LIGHTS)
	    {
	      printf("too many lights, you should increase MAX_LIGHTS!\n");
	      exit(0);
	    }
	  lights[num_lights++] = l;
	}
      else
	{
	  printf("unknown type in scene description:\n%s\n",type);
	  exit(0);
	}
    }
  return 0;
}

//pixel locations are calculated 
void calcPixels(){
	//calculates 4 corner points using equations on slide
	//calculates inner pixel coordinates

	Vertex p1,p2,p3,p4;

    p4.position[0]=(-alpha)*(tan(radFov/2));
    p4.position[1]=tan(radFov/2);
    p4.position[2]=-1.0;
    
    p3.position[0]=(alpha)*(tan(radFov/2));
    p3.position[1]=tan(radFov/2);
    p3.position[2]=-1.0;
    
    p2.position[0]=(alpha)*(tan(radFov/2));
    p2.position[1]=-tan(radFov/2);
    p2.position[2]=-1.0;
    
    p1.position[0]=(-alpha)*(tan(radFov/2));
    p1.position[1]=-tan(radFov/2);
    p1.position[2]=-1.0;

/*the rest of the points are calculated by taking the difference between two corner points, 
 then dividing it by the desired resolution value. using that unit to calculate the position of each
 pixel inbetween*/
	double xUnit = (p2.position[0] - p1.position[0])/resW;
	double yUnit = (p4.position[1] - p1.position[1])/resH;
	  
    for(int j=0;j<resH;j++) 
    {
		pix[0][j].position[0]= p1.position[0]; //starts x index rame for each row
        for(int i=1;i<resW;i++)
        {
            pix[i][j].position[2]=-1.0; //z always the same

			pix[i][j].position[0]= pix[i-1][j].position[0] + xUnit; //x 

			if(j!=0)
			pix[i][j].position[1]= pix[i][j-1].position[1] + yUnit; //y
			else
            pix[i][j].position[1]= p1.position[1]; //y
            
            
        }
    }

}



void display()
{

}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(1,1,1,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
      draw_scene();
      if(mode == MODE_JPEG)
	save_jpg();
    }
  once=1;
}

int main (int argc, char ** argv)
{
  if (argc<2 || argc > 3)
  {  
    printf ("usage: %s <scenefile> [jpegname]\n", argv[0]);
    exit(0);
  }
  //if(argc == 3)
  //  {
      mode = MODE_JPEG;
      filename = "stillImage.jpg";
   //}
  //else if(argc == 2)
   // mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);
  calcPixels(); 
  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}
