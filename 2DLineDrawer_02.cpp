#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;			//The dependences of data structures comes form the standar libraries.

#define xx first		// Redefine the first and second in std::pair<T,T>
#define yy second
#define INF 1e9+5
#define ERR 1e-6			// The margin of aceptable error for work with floating point decimals
typedef long long ll;		
typedef pair<ll,ll> pll;
typedef unsigned char Uchar;
const double PI = acos(-1);			//Constants for rotations
const double RAD = PI/180;


ll xSize = 1920;				// Width of the screen
ll ySize = 1080;				// Height if the screen

double smallizer = 0.5;			// Percentage of scale.


struct triangulo				// This structure is used to save the 3 points in the 2D triangle.
{
	double x[3];
	double y[3];

};

struct line 					// This structure is used for the lines and his processing
{
	double x1, x2, y1, y2;		//The points in the line

	line()						//Default constructor.
	{
		x1 = 0;
		x2 = 0;
		y1 = 0;
		y2 = 0;
	}

};



int setLine(line *ln)		// This funtion evaluate the slope and return 0 if the segment parallel with the y-axis 
{
	if(ln->x1 == ln->x2)
		return 0;
	return 1;
}

int setPoints(line *ln , double x1 , double x2 , double y1 , double y2)	//Discretice a floating points
{
	ln->x1 = x1;
	ln->x2 = x2;
	ln->y1 = y1;
	ln->y2 = y2;
	return setLine(ln);
}

int rotateA(line *ln , double alpha)	// Rotatle alpha angle a segment
{
	double xAux , yAux;

	ln->x1 -= xSize/2;				//Translate the first segment point to the origin
	ln->y1 -= ySize/2;
	xAux = ln->x1;					
	yAux = ln->y1;
	ln->x1 = xAux*cos(alpha) - yAux*sin(alpha);	// Make the rotations
	ln->y1 = yAux*cos(alpha) + xAux*sin(alpha);
	ln->x1 += xSize/2;			// Return the first segment point the original value 
	ln->y1 += ySize/2;


	ln->x2 -= xSize/2;			// The same with the second point
	ln->y2 -= ySize/2;
	xAux = ln->x2;
	yAux = ln->y2;
	ln->x2 = xAux*cos(alpha) - yAux*sin(alpha);
	ln->y2 = yAux*cos(alpha) + xAux*sin(alpha);
	ln->x2 += xSize/2;
	ln->y2 += ySize/2;


	return setLine(ln);			// Return the line rotated
}

bool sortArray(const pll &p1 , const pll &p2)	// The comparation funtion used in std::sort
{
	return p1.xx < p2.xx || (p1.xx == p2.xx && p1.yy < p2.yy);
}


void Bresehand(line *ln , ll flag , vector<pll > &array)		// Get an array of pairs x , y elements in the evaluation for the line
{	
	if(!flag)
	{
		array.assign(abs((ll)(ln->y1-ln->y2)),pll((ll)ln->x1,(ll)min(ln->y1,ln->y2)));
		for(int i = 0 ; i < array.size() ; i++)
			array[i].yy += i;
	}else
	{
		ll x1 = (ll)ln->x1 , x2 = (ll)ln->x2 , y1 = (ll)ln->y1 , y2 = (ll)ln->y2;
		bool negative = false;
		bool pendiente = false;
		if((x1 > x2 && y1 < y2) || (x1 < x2 && y1 > y2))
		{
			swap(y2,y1);
			negative = true;
		}
		ll dx = abs(x1-x2);
		ll dy = abs(y1-y2);
		if(dy == 0)
		{
			for(int i = min(x1,x2) ; i <= max(x1,x2) ; i++)
				array.push_back(pll(i,y1));
			return;
		}
		if(dy > dx)
		{
			pendiente = true;
			swap(x1,y1);
			swap(x2,y2);
			swap(dx,dy);
		}
		if(x2 < x1)
		{
			swap(x1,x2);
			swap(y1,y2);
		}
		bool before = true;	// True is East , False is NorthEast
		ll e = -1;
		for(int i = 0 ; i <= dx ; i++)
		{
			if(e <= 0)
				before = true;
			else
				before = false;
			array.push_back(pll(x1+i,((before)?y1:++y1)));
			if(pendiente)
			{
				swap(array[i].xx,array[i].yy);
			}
			if(before)
			{
				e += 2*dy;
			}
			else
			{
				e += 2*(dy-dx);
			}
			if(i == 0)
			{
				e = 2*dy-dx;
			}
		}
		if(negative)
		{
			for(int i = 0 ; i <= dx/2 ; i++)
			{
				swap(array[i].yy,array[dx-i].yy);
			}
		}
	}

}
vector<pll > makeLineRoast(Uchar r,Uchar g , Uchar b , line *ln , Uchar ***rst , ll flag)	// Put some line in the roaster and return the points that was in the line
{
	vector<pll > array;
	Bresehand(ln,flag,array);
	double aux;

	for(int i = 0 ; i < array.size() ; i++)
	{
		if(array[i].yy >= ySize || array[i].yy < 0 || array[i].xx >= xSize || array[i].xx < 0)
			continue;

		rst[array[i].yy][array[i].xx][0] = r;
		rst[array[i].yy][array[i].xx][1] = g;
		rst[array[i].yy][array[i].xx][2] = b;
	}
	return array;
}

void scanLine(Uchar *** rst ,Uchar r,Uchar g , Uchar b , vector<pll > &array)	// Fill a triangle
{

	for(int i = 0 ; i < array.size()-1 ; i++)
	{
		if(array[i].xx < array[i+1].xx)
			continue;
		else
		{
			for(int j = min(array[i].yy,array[i+1].yy) ; j < max(array[i].yy,array[i+1].yy) ; j++)
			{
				if(array[i].yy >= ySize || array[i].yy < 0 || array[i].xx >= xSize || array[i].xx < 0)	// If the point is outside the screen ignore them
				{
					continue;
				}
				if(rst[j][array[i].xx][0] != (Uchar)255 && rst[j][array[i].xx][1] != (Uchar)255 && rst[j][array[i].xx][2] != (Uchar)255)	// If the pixel is not white ignore them
				{
					continue;
				}
				rst[j][array[i].xx][0] = r;
				rst[j][array[i].xx][1] = g;
				rst[j][array[i].xx][2] = b;
			}
		}
	}
}


void plotRoaster(Uchar *** rst , string fileName)	// Put in a file the result roaster
{

	ofstream ofs;
	ofs.open (fileName, ofstream::out);

	ofs << "P3\n\n"<< xSize << " " << ySize << "\n255\n";
	
	for(int i = ySize-1 ; i >=  0; i--)
	{
		for(int j = 0 ; j < xSize ; j++)
		{
			for(int k = 0 ; k < 3 ; k++)
			{
				ofs << (int)rst[i][j][k] << ' ';
			}
			ofs << " ";
		}
		ofs << "\n";
	}

	ofs.close();
}

void clearRoaster(Uchar *** rst)		// Fill the roaster to white
{
	for(int i = 0 ; i < ySize ; i++)
		for(int j = 0 ; j < xSize ; j++)
			for(int k = 0 ; k < 3 ; k++)
				rst[i][j][k] = 255;
}

int main()
{

	string fName;
	fName = "out.ppm";

	Uchar ***roaster;								// Declare the roaster and assign memory for them
	roaster = new Uchar** [ySize];
	for(int i = 0 ; i < ySize ; i++)
	{
		roaster[i] = new Uchar* [xSize];
		for(int j = 0 ; j < xSize ; j++)
		{
			roaster[i][j] = new Uchar [3];
			for(int k = 0 ; k < 3 ; k++)
				roaster[i][j][k] = 255;
		}
	}
		
	ll flag;
	double minx = 1e9 , maxx = -1e9 , scalex;	// minx save the minimum value for a point in x axis in a all picture
	double miny = 1e9 , maxy = -1e9 , scaley;	// miny is the same but in the y axis
	double maux , basura;						// scalex is a razon to scale the width to the screen size
	vector<triangulo> lineas;					// scaley The same but with the height
	triangulo t1;								// the lineas vector has all the lines in the picture, t1 is an auxiliar triangle


	while(cin >> t1.x[0] >> basura >> t1.y[0] )		// Read the triangles and assing the values for the previous variables
	{
		cin >> t1.x[1] >> basura >> t1.y[1];
		cin >> t1.x[2] >> basura >> t1.y[2];

		if(t1.x[0] > t1.x[1])
			maux = t1.x[1];
		else
			maux = t1.x[0];
		if(maux > t1.x[2])
			maux = t1.x[2];
		if(minx > maux)
			minx = maux;

		if(t1.y[0] > t1.y[1])
			maux = t1.y[1];
		else
			maux = t1.y[0];
		if(maux > t1.y[2])
			maux = t1.y[2];
		if(miny > maux)
			miny = maux;

		maxx = max(max(max(t1.x[0],t1.x[1]),t1.x[2]),maxx);
		maxy = max(max(max(t1.y[0],t1.y[1]),t1.y[2]),maxy);
		
		lineas.push_back(t1);
	}

	if(fabs(maxx - minx) > ERR)							
		scalex = (xSize*smallizer)/(maxx-minx);
	else
		scalex = (xSize * smallizer);
	if(fabs(maxy - miny) > ERR)
		scaley = (ySize * smallizer)/(maxy-miny);
	else
		scaley = (ySize * smallizer);

	for(int i = 0 ; i < lineas.size() ; i++)			// Translate and reescalate the triangles
	{
		for(int j = 0 ; j < 3 ; j++)
		{

			lineas[i].x[j] -= minx;
			lineas[i].y[j] -= miny;
			lineas[i].x[j] *= scalex;
			lineas[i].y[j] *= scaley;
			lineas[i].x[j] += xSize/4;
			lineas[i].y[j] += ySize/4;
		}
	}

	Uchar r,g,b,r2,g2,b2;		//Color variables
	line lnaux;


	r = 0;
	r2 = 198;
	g = 0;
	g2 = 124;
	b = 0;
	b2 = 33;
/*
	for(double k = 0 ; k < 360 ; k += 1.)	// Used only for rotation of the object
	{
		if(k <= 120 && k >= 0)					//
		{										//	This block is used only to change
			r += 2;								//	the color of the line
		}else if(k <= 240 && k > 120)
		{
			r -= 2;
			g += 2;
		}else
		{
			g -= 2;
			b += 2;
		}
*/

		for(int i = 0 ; i < lineas.size() ; i++)
		{
			vector<pll > Triagles;
			vector<pll > addTriague;
			for(int j = 0 ; j < 3 ; j++)
			{
				
				flag = setPoints(&lnaux , lineas[i].x[j] , lineas[i].x[(j+1)%3] , lineas[i].y[j] , lineas[i].y[(j+1)%3]);
		//		if(k > 0)							
		//			flag = rotateA(&lnaux,k*RAD);	
				addTriague = makeLineRoast(r,g,b,&lnaux,roaster,flag);
				Triagles.reserve(addTriague.size() + Triagles.size());
				Triagles.insert(Triagles.end(),addTriague.begin() , addTriague.end());
			}

			sort(Triagles.begin(),Triagles.end(),sortArray);
			scanLine(roaster,r2,g2,b2,Triagles);
		}
		
//		fName = "./OUT/out" + to_string((int)(k + 1.0)) + ".ppm";
		plotRoaster(roaster,fName);
		clearRoaster(roaster);
//	}


	

	for(int i = 0 ; i < ySize ; i++)
	{
		for(int j = 0 ; j < xSize ; j++)
		{
			free(roaster[i][j]);
		}
		free(roaster[i]);
	}
	free(roaster);

	return 0;
}

