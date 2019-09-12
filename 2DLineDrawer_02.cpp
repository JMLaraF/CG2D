#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <map>
#include <set>
#include <queue>

using namespace std;			//The dependences of data structures comes form the standar libraries.

#define xx first		// Redefine the first and second in std::pair<T,T>
#define yy second
#define INF (ll)1<<31
#define FINF 1e9
#define ERR 1e-3			// The margin of aceptable error for work with floating point decimals
typedef long long ll;		
typedef pair<ll,ll> pll;
typedef vector<ll> vll;
typedef vector<pll> vpll;
typedef vector<vll> vvll;
typedef double ld;
typedef vector<ld> vld;
typedef vector<vld> vvld;
typedef unsigned char Uchar;
const double PI = acos(-1);			//Constants for rotations
const double RAD = PI/180;


ld drawERR = 0.015;

double smallizer = 1;			// Percentage of scale.
double cameraDistance = 5;


struct point
{
	double x , y , z;
	point():x(0),y(0),z(0){}
	point(const point &p):x(p.x),y(p.y),z(p.z){}
	point(const double _x , const double _y , const double _z):x(_x),y(_y),z(_z){}
	point operator =(const point &_p)
	{
		this->x = _p.x;
		this->y = _p.y;
		this->z = _p.z;
		return *this;
	}
	point &operator +=(const point &p)
	{
		this->x += p.x;
		this->y += p.y;
		this->z += p.z;
	}
	point &operator -=(const point &p)
	{
		this->x -= p.x;
		this->y -= p.y;
		this->z -= p.z;
	}
	point &operator *=(const double &d)
	{
		this->x *= d;
		this->y *= d;
		this->z *= d;
	}

	point cruz(const point &p)
	{
		double i = this->y*p.z-this->z*p.y;
		double j = this->z*p.x-this->x*p.z;
		double k = this->x*p.y-this->y*p.x;
		return point(i,j,k);
	}
	double punto(const point &p)
	{
		return (this->x*p.x)+(this->y*p.y)+(this->z*p.z);
	}

	void matrizTransform(const vvld &m)
	{
		point aux = (*this);
		this->x = aux.x*m[0][0] + aux.y*m[0][1] + aux.z*m[0][2];
		this->y = aux.x*m[1][0] + aux.y*m[1][1] + aux.z*m[1][2];
		this->z = aux.x*m[2][0] + aux.y*m[2][1] + aux.z*m[2][2];
	}
	double mag()
	{
		return sqrt(this->x*this->x+this->y*this->y+this->z*this->z);
	}
};

point operator +(const point &p1 , const point &p2)
{
	return point(p1)+=p2;
}
point operator -(const point &p1 , const point &p2)
{
	return point(p1)-=p2;
}
point operator *(const point &p , const double &d)
{
	return point(p)*=d;
}
point operator *(const double &d , const point &p)
{
	return point(p)*=d;
}

bool operator ==(const point &p1 , const point &p2)
{
	return (fabs(p1.x - p2.x) <= ERR) && (fabs(p1.y - p2.y) <= ERR) && (fabs(p1.z - p2.z) <= ERR);
}
bool operator !=(const point &p1 , const point &p2)
{
	return (fabs(p1.x - p2.x) > ERR) || (fabs(p1.y - p2.y) > ERR) || (fabs(p1.z - p2.z) > ERR);
}

double cosP(point &p1 , point &p2)
{
	return p1.punto(p2)/(p1.mag()*p2.mag());
}

typedef vector<point> poligono;

struct segment
{
	point p1 , p2;
	segment():p1(),p2(){}
	segment(const point &_p1 , const point &_p2):p1(_p1),p2(_p2){}
	void assign(double x1 , double y1 , double z1 , double x2 , double y2 , double z2)
	{
		p1.x = x1;
		p1.y = y1;
		p1.z = z1;
		p2.x = x2;
		p2.y = y2;
		p2.z = z2;
	}
	void assign(const point &_p1 , const point &_p2)
	{
		p1 = _p1;
		p2 = _p2;
	}
	segment translate(const point &p)
	{
		return segment(this->p1 + p , this->p2 + p);
	}
	bool isVertical()
	{
		return (fabs(p1.x-p2.x) < ERR);
	}
};

//////////////////////////////

struct cara
{
	vll edges;
	point normal;
	point center;
	cara()
	{
		edges.assign(3,-1);
		normal = point(0,0,0);
		center = point(-1,-1,-1);
	}
	cara(const vll &v)
	{
		edges = v;
	}
};




////////////////Roaster/////////////////

struct roaster
{
	ll xSize , ySize;
	vvld zBuffer;
	vector<vvll> rst;
	roaster(ll _x , ll _y)
	{
		xSize = _x;
		ySize = _y;
		zBuffer.assign(_y,vld(_x,FINF));
		rst.assign(_y,vvll(_x,vll(3,255)));
	}
	pll centerOfScreen()
	{
		return pll(xSize/2,ySize/2);
	}
	void clear()
	{
		zBuffer.assign(ySize,vld(xSize,FINF));
		rst.assign(ySize,vvll(xSize,vll(3,255)));
	}
	void resize(ll x, ll y)
	{
		xSize = x;
		ySize = y;
		zBuffer.assign(y,vld(x,FINF));
		rst.assign(y,vvll(x,vll(3,255)));
	}
	bool plot(string fileName)
	{
		ofstream ofs;
		ofs.open (fileName, ofstream::out);
		if(ofs.is_open())
		{
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
			return true;
		}
		return false;
	}
};



/////////////////////////////////////

///////////Proyeccion////////////// Structure maked to conservate the original values of lines and the interpolations of every one of each pixel

struct puntoProyectado
{
	ll x , y , z;
	double ox , oy , oz; //(ox == original x)
	puntoProyectado():x(0),y(0),z(0),ox(0.0),oy(0.0),oz(0.0){}
	puntoProyectado(ll _x , ll _y , ll _z , double _ox , double _oy , double _oz)
	{
		x = _x;
		y = _y;
		z = _z;
		ox = _ox;
		oy = _oy;
		oz = _oz;
	}
	puntoProyectado(const puntoProyectado &p)
	{
		x = p.x;
		y = p.y;
		z = p.z;
		ox = p.ox;
		oy = p.oy;
		oz = p.oz;
	}
};

bool cmpPuntoProyectado(const puntoProyectado &p1 , const puntoProyectado &p2)
{
	return ((p1.x < p2.x)||(p1.x == p2.x && p1.y < p2.y));
}

///////////////////////////////////

/////////////////SpotLigth/////////////////////////////

struct SpotLight
{
	double intensisty;
	double a , b , c , d;
	point pos;
	vector<Uchar> color;
	SpotLight()
	{
		a = 1.0 , b = 1.0 , c = 1.0;
		intensisty = 1.0;
		pos = point(0,0,0);
		color.assign(3,(Uchar)255);
	}
	SpotLight(double x , double y , double z)
	{
		a = 1.0 , b = 1.0 , c = 1.0;
		intensisty = 1.0;
		pos = point(x,y,z);
		color.assign(3,(Uchar)255);
	}
	void setAttenuationConstants(double _a , double _b , double _c)
	{
		a = _a;
		b = _b;
		c = _c;
	}
	void setLightColor(Uchar r , Uchar g , Uchar b)
	{
		color[0] = r;
		color[1] = g;
		color[2] = b;
	}
	void setLigthPosition(double _x , double _y , double _z)
	{
		pos = point(_x,_y,_z);
	}
};

//////////////////////////////////////////////////////

double pot(double a , ll b)
{
	double ans = 1.0;
	for(int i = 0 ; i < 60 ; i++)
	{
		if((((ll)1<<i)&b)!=0)
			ans *= a;
		a *= a;
	}
	return ans;
}

vector<Uchar> setLigth(const SpotLight &spl ,cara &Cara ,const vector<Uchar> color)
{
	vector<Uchar> ans(3,0);
	double Lr , Lg , Lb , d , cosAngle;
	point pAux = spl.pos-Cara.center;
	d = pAux.mag();
	cosAngle = max(cosP(Cara.normal,pAux),0.0);
	double attenuation = 1.0/(spl.a*d*d + spl.b*d + spl.c);
	for(int i = 0 ; i < ans.size() ; i++)
		ans[i] += (Uchar)min(255.0,((double)color[i]+((double)spl.color[i]*attenuation*cosAngle)));
	return ans;
}

vll specularLigthing(const point &pixel , const SpotLight &spl ,const point &normal ,point &camVector)
{
	point d = pixel - spl.pos;
	point r = d - normal*(2*(d.punto(normal)));
	double dis = r.mag();
	double attenuation = 1.0/(spl.a*dis*dis + spl.b*dis + spl.c);
	double cosAngle = cosP(r,camVector);
//	cosAngle *= -1.0;
	cosAngle = max(cosAngle*(-1),0.0);
	cosAngle = pow(cosAngle,3.9);
	double factor = cosAngle;
	vll ans ={(ll)(((double)spl.color[0])*factor),(ll)(((double)spl.color[1])*factor),(ll)(((double)spl.color[2])*factor)};
	return ans;

}


bool sortArray(const pair<pll,ll> &p1 , const pair<pll,ll> &p2)	// The comparation funtion used in std::sort
{
	return p1.xx.xx < p2.xx.xx || (p1.xx.xx == p2.xx.xx && p1.xx.yy < p2.xx.yy);
}


void Bresehand(segment &ln , vector<puntoProyectado> &array , segment &originalSegment)		// Get an array of pairs x , y elements in the evaluation for the line
{	
	ll x1 = (ll)ln.p1.x , x2 = (ll)ln.p2.x , y1 = (ll)ln.p1.y , y2 = (ll)ln.p2.y , z1 = ln.p1.z , z2 = ln.p2.z;
	ll minZ = min(z1,z2);
	ll dx = abs(x2-x1) , dy = abs(y2-y1) , dz = abs(z2-z1);
	double fdx = fabs(originalSegment.p2.x-originalSegment.p1.x) , fdy = fabs(originalSegment.p2.y-originalSegment.p1.y) , fdz = fabs(originalSegment.p2.z-originalSegment.p1.z);
	if(ln.isVertical())
	{
		int i , j , k = 0 , dj;
		double I , J , dI , dJ;
		
		if(y1 < y2)
		{
			i = y1;
			I = originalSegment.p1.y;
			j = z1;
			J = originalSegment.p1.z;
			
		}else
		{
			i = y2;
			I = originalSegment.p2.y;
			j = z2;
			J = originalSegment.p2.z;
		}
		if(dy == 0)
		{
			array.push_back(puntoProyectado(x1,i,j,originalSegment.p1.x,originalSegment.p1.y,originalSegment.p1.z));
			return;
		}
		for(; i <= max(y1,y2) ; i++ , k++)
		{
			dj = (dz/dy)*k;
			dI = (fdy/dy)*k;
			dJ = (fdz/dy)*k;
			array.push_back(puntoProyectado(x1,i,j+((j == min(z1,z2)?dj:-dj)),originalSegment.p1.x,I+((I == min(originalSegment.p1.y,originalSegment.p2.y))?dI:-dI),J+((J == min(originalSegment.p1.z,originalSegment.p2.z))?dJ:-dJ)));
		}
	}else
	{
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
			if(dx == 0)
			{
				array.push_back(puntoProyectado(x1,y2,z1,originalSegment.p1.x,originalSegment.p1.y,originalSegment.p1.z));
				return;
			}
			int i , j , k = 0 , dj;
			double I , J , dI , dJ;
			
			if(x1 < x2)
			{
				i = x1;
				I = originalSegment.p1.x;
				j = z1;
				J = originalSegment.p1.z;
				
			}else
			{
				i = x2;
				I = originalSegment.p2.x;
				j = z2;
				J = originalSegment.p2.z;
			}
			for(; i <= max(x1,x2) ; i++ , k++)
			{
				dj = (dz/dx)*k;
				dI = (fdx/dx)*k;
				dJ = (fdz/dx)*k;
				array.push_back(puntoProyectado(i,y1,j+((j == min(z1,z2)?dj:-dj)),I+((I == min(originalSegment.p1.x,originalSegment.p2.x))?dI:-dI),originalSegment.p1.y,J+((J == min(originalSegment.p1.z,originalSegment.p2.z))?dJ:-dJ)));
			}
			return;
		}
		
		if(dy > dx)
		{
			pendiente = true;
			swap(x1,y1);
			swap(x2,y2);
			swap(dx,dy);
			swap(originalSegment.p1.x,originalSegment.p1.y);
			swap(originalSegment.p2.x,originalSegment.p2.y);
			swap(fdx,fdy);
		}

		if(x2 < x1)
		{
			swap(x1,x2);
			swap(y1,y2);
			swap(originalSegment.p1.x,originalSegment.p2.x);
			swap(originalSegment.p1.y,originalSegment.p2.y);
		}
		double I = originalSegment.p1.x , J = originalSegment.p1.z , K = originalSegment.p1.y;
		ll j = z1;
		bool before = true;	// True is East , False is NorthEast
		ll e = -1;
		for(int i = 0 ; i <= dx ; i++)
		{
			ll dj = (dz/dx)*i;
			double dI = (fdx/dx)*i , dJ = (fdz/dx)*i , dK = (fdy/dx)*i;
			if(e <= 0)
				before = true;
			else
				before = false;
			array.push_back(puntoProyectado(x1+i,((before)?y1:++y1),minZ+(ll)(round(fdz/dx)*i),I+((I == min(originalSegment.p1.x,originalSegment.p2.x))?dI:-dI),K+((K == min(originalSegment.p1.y,originalSegment.p2.y))?dK:-dK),J+((J == min(originalSegment.p1.z,originalSegment.p2.z))?dJ:-dJ)));
			if(pendiente)
			{
				swap(array[i].x,array[i].y);
				swap(array[i].ox,array[i].oy);
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
				swap(array[i].y,array[dx-i].y);
				swap(array[i].oy,array[dx-i].oy);
			}
		}
	}
}


vector<puntoProyectado> makeLineRoast(vector<Uchar> color , segment &ln , segment originalSegment , roaster &Roaster)	// Put some line in the roaster and return the points that was in the line
{
	vector<puntoProyectado> array;
	Bresehand(ln,array,originalSegment);
	double aux;
	pll center = Roaster.centerOfScreen();
	for(int i = 0 ; i < array.size() ; i++)
	{

		ll xAux = array[i].x+center.xx , yAux = array[i].y+center.yy;
		if(yAux >= Roaster.ySize || yAux < 0 || xAux >= Roaster.xSize || xAux < 0)
			continue;
		if(array[i].oz >= Roaster.zBuffer[yAux][xAux])
			continue;
	//	Roaster.zBuffer[yAux][xAux] = array[i].oz;
	//	for(int j = 0 ; j < 3 ; j++)
	//		Roaster.rst[yAux][xAux][j] = color[j];
	}
	return array;
}

void scanLine(roaster &Roaster ,vector<Uchar> color, vector<puntoProyectado> &array , SpotLight &spl , point &camVector , point &faceNormal)	// Fill a triangle
{
	if(array.empty())
		return;
	sort(array.begin(),array.end(),cmpPuntoProyectado);
	pll center = Roaster.centerOfScreen();
	for(int i = 0 ; i < array.size()-1 ; i++)
	{
		if(array[i].x < array[i+1].x)
			continue;
		else
		{
			ll xAux = array[i].x+center.xx;
			double fdx = fabs(array[i].ox-array[i+1].ox), dy = abs(array[i].y-array[i+1].y) , fdy = fabs(array[i].oy-array[i+1].oy) , fdz = fabs(array[i].oz-array[i+1].oz);
			double xIndex , yIndex , zIndex;
			bool dir = false;
			if(array[i+1].y == min(array[i].y,array[i+1].y))
				dir = true;
			for(int j = min(array[i].y,array[i+1].y)+center.yy , k = 0 ; j <= max(array[i].y,array[i+1].y)+center.yy ; j++ , k++)
			{
				 
				if(j >= Roaster.ySize ||j < 0 || xAux >= Roaster.xSize || xAux < 0)	// If the point is outside the screen ignore them
					continue;
			//	if(Roaster.rst[j][xAux][0] != (Uchar)255 && Roaster.rst[j][xAux][1] != (Uchar)255 && Roaster.rst[j][xAux][2] != (Uchar)255)	// If the pixel is not white ignore them
			//		continue;
			
				if(dir)
				{
					xIndex = array[i+1].ox + ((array[i+1].ox == min(array[i].ox,array[i+1].ox))?(fdx/dy)*k:-(fdx/dy)*((double)k));
					yIndex = array[i+1].oy + ((array[i+1].oy == min(array[i].oy,array[i+1].oy))?(fdy/dy)*k:-(fdy/dy)*((double)k));
					zIndex = array[i+1].oz + ((array[i+1].oz == min(array[i].oz,array[i+1].oz))?(fdz/dy)*k:-(fdz/dy)*((double)k));
				}else
				{
					xIndex = array[i].ox + ((array[i].ox == min(array[i].ox,array[i+1].ox))?(fdx/dy)*k:-(fdx/dy)*(k));
					yIndex = array[i].oy + ((array[i].oy == min(array[i].oy,array[i+1].oy))?(fdy/(double)dy)*(double)k:-(fdy/(double)dy)*((double)k));
					zIndex = array[i].oz + ((array[i].oz == min(array[i].oz,array[i+1].oz))?(fdz/(double)dy)*(double)k:-(fdz/(double)dy)*((double)k));
				}
				if(zIndex >= Roaster.zBuffer[j][xAux])
					continue;
				Roaster.zBuffer[j][xAux] = zIndex;
				vll specularLigth = specularLigthing(point(xIndex,yIndex,zIndex),spl,faceNormal,camVector);
				for(int k = 0 ; k < 3 ; k++)
					Roaster.rst[j][xAux][k] = (Uchar)min((ll)254,(((ll)color[k])+specularLigth[k]));
			}
		}
	}
}



////////////////////Main funtions ///////////////////////



//////////////////////////////////////////////////////////



/////////////////VLF Format /////////////////////////////

class Object
{
	private:
		void readRawFormat(fstream &input)
		{
			map<pair<pair<ld,ld>,ld>,ll> listOfVectors;
			map<pll,ll> listOfEdges;
			pair<pair<ld,ld>,ld> vAux[3];
			while(input >> vAux[0].xx.xx >> vAux[0].xx.yy >> vAux[0].yy)
			{
				input >> vAux[1].xx.xx >> vAux[1].xx.yy >> vAux[1].yy;
				input >> vAux[2].xx.xx >> vAux[2].xx.yy >> vAux[2].yy;
				for(int i = 0 ; i < 3 ; i++)
				{
					if(!listOfVectors.count(vAux[i]))
					{
						listOfVectors[vAux[i]] = vectores.size();
						vectores.push_back(point(vAux[i].xx.xx,vAux[i].xx.yy,vAux[i].yy));
					}
				}
				ll a = listOfVectors[vAux[0]] , b = listOfVectors[vAux[1]] , c = listOfVectors[vAux[2]];
				ll x , y , z;
				if(!listOfEdges.count(pll(a,b)) && !listOfEdges.count(pll(b,a)))
				{
					listOfEdges[pll(a,b)] = aristas.size();
					x = aristas.size();
					aristas.push_back(pll(a,b));
				}else
				{
					x = (listOfEdges.count(pll(a,b)))?listOfEdges[pll(a,b)]:listOfEdges[pll(b,a)];
				}

				if(!listOfEdges.count(pll(b,c)) && !listOfEdges.count(pll(c,b)))
				{
					listOfEdges[pll(b,c)] = aristas.size();
					y = aristas.size();
					aristas.push_back(pll(b,c));
				}else
				{
					y = (listOfEdges.count(pll(b,c)))?listOfEdges[pll(b,c)]:listOfEdges[pll(c,b)];
				}

				if(!listOfEdges.count(pll(c,a)) && !listOfEdges.count(pll(a,c)))
				{
					listOfEdges[pll(c,a)] = aristas.size();
					z = aristas.size();
					aristas.push_back(pll(c,a));
				}else
				{
					z = (listOfEdges.count(pll(c,a)))?listOfEdges[pll(c,a)]:listOfEdges[pll(a,c)];
				}
				caras.push_back(cara(vll(3)));
				caras[caras.size()-1].edges[0] = x;
				caras[caras.size()-1].edges[1] = y;
				caras[caras.size()-1].edges[2] = z;
			}
		}
		void setNormals()
		{
			pll paux;
			point vecAux[3];
			ll aux;
			vpll direcciones(aristas.size(),pll(0,0));
			vll e;
			set<ll> carasUsadas;
			queue<ll> cola;
			cola.push(0);
			carasUsadas.insert(0);
		//	direcciones[caras[0].edges[0]].xx = 1;
			while(!cola.empty())
			{
				aux = cola.front();
				cola.pop();
				e = caras[aux].edges;
				for(int i = 0 ; i < caras.size() ; i++)
				{
					for(int j = 0 ; j < 3 ; j++)
					{
						if(!carasUsadas.count(i) && (caras[i].edges[j] == e[0] || caras[i].edges[j] == e[1] || caras[i].edges[j] == e[2]))
						{
							cola.push(i);
							carasUsadas.insert(i);
						}
					}
				}
				ll k , h;
				bool isFirstUsed = true;
				for(int i = 0 ; i < 3 ; i++)
				{
					if(direcciones[e[i]].xx != 0)
					{
						direcciones[e[i]].yy = (-1)*direcciones[e[i]].xx;
						h = direcciones[e[i]].yy;
						paux = ((h > 0)? aristas[e[i]] : pll(aristas[e[i]].yy,aristas[e[i]].xx));
						k = i;
						isFirstUsed = false;
						break;
					}
				}
				if(isFirstUsed)
				{
					direcciones[e[0]].xx = 1;
					h = direcciones[e[0]].xx;
					paux = ((h > 0)? aristas[e[0]] : pll(aristas[e[0]].yy,aristas[e[0]].xx));
					k = 0;
				}
				vecAux[0] = vectores[paux.xx];
				vecAux[1] = vectores[paux.yy];
				if(aristas[e[(k+1)%3]].xx == paux.yy)
				{
					vecAux[2] = vectores[aristas[e[(k+1)%3]].yy];
					direcciones[e[(k+1)%3]].xx = 1;
					if(aristas[e[(k+2)%3]].xx == aristas[e[(k+1)%3]].yy)
						direcciones[e[(k+2)%3]].xx = 1;
					else
						direcciones[e[(k+2)%3]].xx = -1;
				}
				else if(aristas[e[(k+2)%3]].xx == paux.yy)
				{
					vecAux[2] = vectores[aristas[e[(k+2)%3]].yy];
					direcciones[e[(k+2)%3]].xx = 1;
					if(aristas[e[(k+1)%3]].xx == aristas[e[(k+2)%3]].yy)
						direcciones[e[(k+1)%3]].xx = 1;
					else
						direcciones[e[(k+1)%3]].xx = -1;
				}else if(aristas[e[(k+1)%3]].yy == paux.yy)
				{
					vecAux[2] = vectores[aristas[e[(k+1)%3]].xx];
					direcciones[e[(k+1)%3]].xx = -1;
					if(aristas[e[(k+2)%3]].xx == aristas[e[(k+1)%3]].xx)
						direcciones[e[(k+2)%3]].xx = 1;
					else
						direcciones[e[(k+2)%3]].xx = -1;
				}
				else if(aristas[e[(k+2)%3]].yy == paux.yy)
				{
					vecAux[2] = vectores[aristas[e[(k+2)%3]].xx];
					direcciones[e[(k+2)%3]].xx = -1;
					if(aristas[e[(k+1)%3]].xx == aristas[e[(k+2)%3]].xx)
						direcciones[e[(k+1)%3]].xx = 1;
					else
						direcciones[e[(k+1)%3]].xx = -1;
				}
				point p1 = vecAux[0]-vecAux[1];
				point p2 = vecAux[2]-vecAux[1];
				caras[aux].normal = p1.cruz(p2);
				caras[aux].normal *= (1/caras[aux].normal.mag());
				ll cond = cola.size();
				if(cola.empty())
				{
					for(int i = 0 ; i < caras.size() ; i++)
					{
						if(!carasUsadas.count(i))
						{
							cola.push(i);
							carasUsadas.insert(i);
							direcciones[caras[i].edges[0]].xx = 1;
							break;
						}
					}
				}
			}
		}
		void setCenters()
		{
			for(int i = 0 ; i < caras.size() ; i++)
			{
				point pAux[3];
				pAux[0] = vectores[aristas[caras[i].edges[0]].xx];
				pAux[1] = vectores[aristas[caras[i].edges[0]].yy];
				if(vectores[aristas[caras[i].edges[1]].xx] != pAux[0] && vectores[aristas[caras[i].edges[1]].xx] != pAux[1])
					pAux[2] = vectores[aristas[caras[i].edges[1]].xx];
				else
					pAux[2] = vectores[aristas[caras[i].edges[1]].yy];
				caras[i].center = point((pAux[0].x+pAux[1].x+pAux[2].x)/3,(pAux[0].y+pAux[1].y+pAux[2].y)/3,(pAux[0].z+pAux[1].z+pAux[2].z)/3);
			}
		}

	public:
		vector<point> vectores;
		vpll aristas;
		vector<cara> caras;
		Object(const string &s , char format)
		{
			fstream input;
			input.open(s);
			if(input.is_open())
			{
				switch(format)
				{
					case 'R':
						readRawFormat(input);
						setNormals();
						setCenters();
						for(int i = 0 ; i < caras.size() ; i++)
							caras[i].normal *= (-1);
					break;
				}
			}else
				cout << "Error al abrir archivo\n";
			input.close();
		}

};


/////////////////////////////////////////////////////////


///////////////// Camera ///////////////////////////////

class camera
{
	private:
		double d;
		point pos;
		point upVec;
		point dirVec;
		vvld matrizDeRotacion;
		vvld matrizDeProyeccion;

		point proyectar(const point &p)
		{
			double a = p.x*matrizDeProyeccion[0][0] + p.y*matrizDeProyeccion[0][1] + p.z*matrizDeProyeccion[0][2] + matrizDeProyeccion[0][3];
			double b = p.x*matrizDeProyeccion[1][0] + p.y*matrizDeProyeccion[1][1] + p.z*matrizDeProyeccion[1][2] + matrizDeProyeccion[1][3];
			double c = p.x*matrizDeProyeccion[2][0] + p.y*matrizDeProyeccion[2][1] + p.z*matrizDeProyeccion[2][2] + matrizDeProyeccion[2][3];
			double w = p.x*matrizDeProyeccion[3][0] + p.y*matrizDeProyeccion[3][1] + p.z*matrizDeProyeccion[3][2] + matrizDeProyeccion[3][3];
			if(fabs(w) <= ERR)
				w = 0.0001;
			return point(a/w,b/w,c/w);
		}
	public:
		camera()
		{
			pos = point(0,0,0);
			dirVec = point(0,0,1);
			d = 2030;
			matrizDeRotacion.assign(3,vld(3,0));
			matrizDeProyeccion.assign(4,vld(4,0));
			matrizDeProyeccion[0][0] = 1;
			matrizDeProyeccion[1][1] = 1;
			matrizDeProyeccion[2][2] = 1;
			matrizDeProyeccion[3][2] = 1.0/d;
			matrizDeRotacion[0][0] = 1;
			matrizDeRotacion[1][1] = 1;
			matrizDeRotacion[2][2] = 1;
		}
		camera(double a , double b , double c , double x , double y , double z)
		{
			pos = point(x,y,z);
			dirVec = point(0,0,1);
			d = 5;
			matrizDeRotacion.assign(3,vld(3,0));
			matrizDeProyeccion.assign(4,vld(4,0));
			matrizDeProyeccion[0][0] = 1;
			matrizDeProyeccion[1][1] = 1;
			matrizDeProyeccion[2][2] = 1;
			matrizDeProyeccion[3][2] = 1.0/d;
			matrizDeRotacion[0][0] = 1;
			matrizDeRotacion[1][1] = 1;
			matrizDeRotacion[2][2] = 1;
			rotateDir(a,b,c);
		}
		void setDistance(double _d)
		{
			d = _d;
			matrizDeProyeccion[3][2] = 1.0/d;
		}
		void Translate(double _x , double _y , double _z)
		{
			pos.x = _x;
			pos.y = _y;
			pos.z = _z;
		}
		void rotateDir (double a , double b , double c)
		{
			double ca = cos(a) , cb = cos(b) , cc = cos(c) , sa = sin(-a) , sb = sin(-b) , sc = sin(-c);
			if(a == 0)
				sa = 0;
			if(b == 0)
				sb = 0;
			if(c == 0)
				sc = 0;
			matrizDeRotacion[0][0] = cb*cc;
			matrizDeRotacion[0][1] = ((sc==0)?0:-cb*sc);
			matrizDeRotacion[0][2] = sb;
			matrizDeRotacion[1][0] = sa*sb*cc+ca*sc;
			matrizDeRotacion[1][1] = ca*cc-(sa*sb*sc);
			matrizDeRotacion[1][2] = ((sa==0)?0:-sa*cb);
			matrizDeRotacion[2][0] = sa*sc-ca*cc*sb;
			matrizDeRotacion[2][1] = sa*cc+sb*ca*sc;
			matrizDeRotacion[2][2] = ca*cb;

			vvld mz = matrizDeRotacion;
			sa = sin(a) , sb = sin(b) , sc = sin(c);
			mz[0][0] = cb*cc;
			mz[0][1] = ((sc==0)?0:-cb*sc);
			mz[0][2] = sb;
			mz[1][0] = sa*sb*cc+ca*sc;
			mz[1][1] = ca*cc-(sa*sb*sc);
			mz[1][2] = ((sa==0)?0:-sa*cb);
			mz[2][0] = sa*sc-ca*cc*sb;
			mz[2][1] = sa*cc+sb*ca*sc;
			mz[2][2] = ca*cb;

			dirVec = point(0,0,1);
			dirVec.matrizTransform(mz);
//			cout << dirVec.x << ' ' << dirVec.y << ' ' << dirVec.z << '\n';
			/* 
				display matrizDeRotacion
			*/
		}
		void rotateDir (double a , double b , double c , double _x , double _y , double _z)
		{
			Translate(_x,_y,_z);
			dirVec = dirVec - pos;
			rotateDir(a,b,c);

		}

		void TakeSnap(Object &obj , SpotLight &spl , string fName , roaster &Roaster , vector<Uchar> &lineColor , vector<Uchar> &fillColor)
		{
			ll aux;
			segment lnaux , lnOrigin;
			bool flag;
			point pointAux[4];
			vector<puntoProyectado> array;
			vector<puntoProyectado> filled;
			for(int i = 0 ; i < obj.caras.size() ; i++)
			{
		//		cout << "Nor = " << obj.caras[i].normal.x << ' ' << obj.caras[i].normal.y << ' ' << obj.caras[i].normal.z << " | " << cosP(dirVec , obj.caras[i].normal) << '\n';
		//		if(obj.caras[i].normal.z > 0.5 && obj.caras[i].normal.z <= 1.0)
		//			cout << i << " <I\n";
				if(cosP(dirVec , obj.caras[i].normal) > 1.8 )
					continue;
	
				for(int j = 0 ; j < 3 ; j++)
				{
					aux = obj.caras[i].edges[j];
					pointAux[0] = obj.vectores[obj.aristas[aux].xx];
					pointAux[1] = obj.vectores[obj.aristas[aux].yy];
					pointAux[0] -= this->pos;
					pointAux[1] -= this->pos;
					pointAux[0].matrizTransform(this->matrizDeRotacion);
					pointAux[1].matrizTransform(this->matrizDeRotacion);
					if(pointAux[0].z < 1.0 && pointAux[1].z < 1.0)
						continue;
					pointAux[2] = proyectar(pointAux[0]);
					pointAux[3] = proyectar(pointAux[1]);
					pointAux[2].z = pointAux[0].z;
					pointAux[3].z = pointAux[1].z;
					lnaux.assign(pointAux[2],pointAux[3]);
					lnOrigin.assign(obj.vectores[obj.aristas[aux].xx],obj.vectores[obj.aristas[aux].yy]);
					array = makeLineRoast(lineColor,lnaux,lnOrigin,Roaster);
					for(int k = 0 ; k < array.size() ; k++)
					filled.push_back(array[k]);
				}
				vector<Uchar> LfillColor = setLigth(spl,obj.caras[i],fillColor);
				scanLine(Roaster,LfillColor,filled,spl,this->dirVec,obj.caras[i].normal);
				filled.clear();
			}
			Roaster.plot(fName);
			Roaster.clear();

		}
};

////////////////////////////////////////////////////////

int main()	//Compile with c++11	// g++ 2DLineDrawer_02.cpp -std=c++11
{

	string fName;
	fName = "out.ppm";
	roaster Roaster(1920,1080); 
	vector<Uchar> lineColor(3,(Uchar)0);
	vector<Uchar> fillColor(3,(Uchar)50);

	string s;
	cin >> s;
	Object cubes(s,'R');
	camera cam;
	SpotLight ligth(0.0,6.0,-2.0);
	ligth.setAttenuationConstants(0.0005,0.2,0.0);
	ligth.setLightColor(170,170,0);
	int i = 0;
	for(i = 0 ; i < 40 ; i++)
	{
		s = "./OUT/out0"+to_string(i)+".ppm";
		cout << s << '\n';
		ligth.setLigthPosition(0.0,6.0-((6.0/20.0)*i),2.0);

//		cam.rotateDir(1.0*i*RAD,0.5*i*RAD,0.0);
		cam.TakeSnap(cubes,ligth,s,Roaster,lineColor,fillColor);	
//		ligth.setLightColor(255,255,255);
//		string s2 = "./OUT/out0"+to_string(255)+".ppm";
//		cam.TakeSnap(cubes,ligth,s2,Roaster,lineColor,fillColor);
	}
	



	return 0;
}

