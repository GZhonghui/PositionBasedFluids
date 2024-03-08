#include<GL/glut.h>
#include<GL/glu.h>
#include<GL/gl.h>

#include<algorithm>
#include<vector>
#include<cmath>

#include<cstdio>

using namespace std;

class Vector3;
struct Particle;

const double Pi(acos(-1.0));
const double g (9.8);
const double dt(.05);
const double rh(0.5);
const double h (2.0);
const double r0(100);

const int maxIter(5);
const int screenWidth(1000);
const int screenHeight(800);

const int boxWidth =32;
const int fluidWidth=10;

class Vector3
{
public:
    double x, y, z;

    Vector3() : x(0), y(0), z(0) {}
    Vector3(double xx) : x(xx), y(xx), z(xx) {}
    Vector3(double xx, double yy, double zz) : x(xx), y(yy), z(zz) {}

    Vector3 operator * (const double &r) const { return Vector3(x * r, y * r, z * r); }
    Vector3 operator / (const double &r) const { return Vector3(x / r, y / r, z / r); }

    double norm() const {return sqrt(x * x + y * y + z * z);}
    Vector3 normalized()
    {
        double n = std::sqrt(x * x + y * y + z * z);
        return Vector3(x / n, y / n, z / n);
    }

    Vector3 operator * (const Vector3 &v) const { return Vector3(x * v.x, y * v.y, z * v.z); }
    Vector3 operator - (const Vector3 &v) const { return Vector3(x - v.x, y - v.y, z - v.z); }
    Vector3 operator + (const Vector3 &v) const { return Vector3(x + v.x, y + v.y, z + v.z); }
    Vector3 operator - () const { return Vector3(-x, -y, -z); }

    Vector3& operator += (const Vector3 &v)
    {
        x += v.x, y += v.y, z += v.z;
        return *this;
    }

    friend Vector3 operator * (const float &r, const Vector3 &v)
    {
        return Vector3(v.x * r, v.y * r, v.z * r);
    }

    static Vector3 Min(const Vector3 &p1, const Vector3 &p2)
    {
        return Vector3(std::min(p1.x, p2.x), std::min(p1.y, p2.y),std::min(p1.z, p2.z));
    }
    static Vector3 Max(const Vector3 &p1, const Vector3 &p2)
    {
        return Vector3(std::max(p1.x, p2.x), std::max(p1.y, p2.y),std::max(p1.z, p2.z));
    }
};

struct Particle
{
public:
	bool m_fluids;
	double m_mass;
	double m_lambda;
	double m_density;
	Vector3 m_deltaPos;
	Vector3 m_position;
	Vector3 m_lastPosition;
	Vector3 m_velocity;
	vector<unsigned int> m_neighborhood;
};

vector<Particle> Particles;

double poly6WKernel(const Vector3 & r)
{
	double ret = 0.0;
	double rl = r.norm();
	double q = rl / h;
	double h3 = h*h*h;
	if (q <= 0.5)
	{
		double q2 = q * q;
		double q3 = q2 * q;
		ret = 8.0 / (Pi * h3) * (6.0 * q3 - 6.0 * q2 + 1.0);
	}
	else
	{
		ret = 16.0 / (Pi * h3) * pow(1 - q, 3.0);
	}
	return ret;
}

Vector3 spikyWKernelGrad(const Vector3 & r)
{
	Vector3 ret(0.0);
	double rl = r.norm();
	double q = rl / h;
	double h3 = h*h*h;
	if (rl > 1.0e-6)
	{
		const Vector3 gradq=(1.0 / (rl * h)) * r;
		if (q <= 0.5)
		{
			ret = (48.0 / (Pi * h3) * q * (3.0 * q - 2.0)) * gradq;
		}
		else
		{
			double factor = 1.0 - q;
			ret = (48.0 / (Pi * h3) * (-factor * factor)) * gradq;
		}
	}
	return ret;
}

void updateNeighborhood(bool fluids)
{
	if(fluids)
	{
		for(unsigned int x=0;x<Particles.size();x+=1)
		{
			if(!Particles[x].m_fluids)
				break;
			Particles[x].m_neighborhood.clear();
			for(unsigned y=0;y<Particles.size();y+=1)
			{
				if((Particles[x].m_position-Particles[y].m_position).norm()<h)
				Particles[x].m_neighborhood.push_back(y);
			}
		}
	}else
	{
		for(unsigned int x=Particles.size()-1;x>=0;x-=1)
		{
			if(Particles[x].m_fluids)
				break;
			Particles[x].m_neighborhood.clear();
			for(unsigned y=0;y<Particles.size();y+=1)
			{
				if(!Particles[y].m_fluids&&(Particles[x].m_position-Particles[y].m_position).norm()<h)
				{
					Particles[x].m_neighborhood.push_back(y);
				}
			}
		}
	}
}

void computeFluidDensity(unsigned int particleIndex)
{
	double density = Particles[particleIndex].m_mass * poly6WKernel(Vector3(0));
	for (unsigned int x = 0; x < Particles[particleIndex].m_neighborhood.size(); x+=1)
	{
		unsigned int neighborIndex = Particles[particleIndex].m_neighborhood[x];
		density += Particles[neighborIndex].m_mass * poly6WKernel(Particles[particleIndex].m_position - Particles[neighborIndex].m_position);
	}
	Particles[particleIndex].m_density=density;
}

void computeLagrangeMultiplier(unsigned int particleIndex)
{
	double lambda = 0.0;
	const double eps = 1.0e-6;
	if (Particles[particleIndex].m_density > r0)
	{
		double sum_grad_cj = 0.0;
		Vector3 grad_ci(0.0);
		for (unsigned int x = 0; x < Particles[particleIndex].m_neighborhood.size(); ++x)
		{
			unsigned int neighborIndex = Particles[particleIndex].m_neighborhood[x];
			Vector3 grad_cj = (Particles[neighborIndex].m_mass / r0)
				* spikyWKernelGrad(Particles[particleIndex].m_position - Particles[neighborIndex].m_position);
			sum_grad_cj += pow(grad_cj.norm(), 2.0);
			grad_ci += grad_cj;
		}
		sum_grad_cj += pow(grad_ci.norm(), 2.0);
		lambda = -(Particles[particleIndex].m_density/r0-1.0) / (sum_grad_cj + eps);
	}
	Particles[particleIndex].m_lambda=lambda;
}

void solveDensityConstraint(unsigned int particleIndex)
{
	Vector3 deltaPos = Vector3(0.0f);
	for (unsigned int x = 0; x < Particles[particleIndex].m_neighborhood.size(); ++x)
	{
		unsigned int neighborIndex = Particles[particleIndex].m_neighborhood[x];
		if (Particles[neighborIndex].m_fluids)
		{
			Vector3 grad_cj = (Particles[neighborIndex].m_mass / r0) 
				* spikyWKernelGrad(Particles[particleIndex].m_position - Particles[neighborIndex].m_position);
			deltaPos += (Particles[particleIndex].m_lambda + Particles[neighborIndex].m_lambda) * grad_cj;
		}
		else
		{
			Vector3 grad_cj = (Particles[neighborIndex].m_mass / r0)
				* spikyWKernelGrad(Particles[particleIndex].m_position - Particles[neighborIndex].m_position);
			deltaPos += (Particles[particleIndex].m_lambda) * grad_cj;
		}
	}
	Particles[particleIndex].m_deltaPos=deltaPos;
}

void velocityUpdateFirstOrder(unsigned int particleIndex)
{
	Particles[particleIndex].m_velocity = (1.0 / dt) * (Particles[particleIndex].m_position - Particles[particleIndex].m_lastPosition);
}

void constraintProjection()
{
	unsigned int iter = 0;

	while (iter < maxIter)
	{
		updateNeighborhood(true);

		for (unsigned int x = 0; x < Particles.size(); x+=1)
		{
			if(!Particles[x].m_fluids)
				break;
			computeFluidDensity(x);
			computeLagrangeMultiplier(x);
		}
		
		for (unsigned int x = 0; x < Particles.size(); x+=1)
		{
			if(!Particles[x].m_fluids)
				break;
			solveDensityConstraint(x);
		}

		for (unsigned int x = 0; x < Particles.size(); x+=1)
		{
			if(!Particles[x].m_fluids)
				break;
			Particles[x].m_position += Particles[x].m_deltaPos;
			Particles[x].m_position.x=max(1.,min(boxWidth-1.,Particles[x].m_position.x));
			Particles[x].m_position.y=max(1.,min(boxWidth-1.,Particles[x].m_position.y));
			Particles[x].m_position.z=max(1.,min(boxWidth-1.,Particles[x].m_position.z));
		}
		++iter;
	}
}

void applyGravity(unsigned int particleIndex)
{
	Vector3 gravity(0,-g,0);
	Particles[particleIndex].m_velocity+=gravity*dt;
	Particles[particleIndex].m_position+=Particles[particleIndex].m_velocity*dt;
}

void computeXSPHViscosity()
{
	double c=0.02;
	for (unsigned int x = 0; x < Particles.size(); x+=1)
	{
		Vector3 velocity = Particles[x].m_velocity;
		Vector3 position = Particles[x].m_position;
		Vector3 sum_value(0.0);
		for (unsigned int y = 0; y < Particles[x].m_neighborhood.size(); y+=1)
		{
			unsigned int neighborIndex = Particles[x].m_neighborhood[y];
			if(Particles[neighborIndex].m_fluids)
			{
				const double density_j = Particles[neighborIndex].m_density;
				const Vector3 position_j = Particles[neighborIndex].m_position;
				const Vector3 velocity_j = Particles[neighborIndex].m_velocity;
				Vector3 tmp = velocity - velocity_j;
				tmp = tmp * poly6WKernel(position - position_j) * (Particles[neighborIndex].m_mass / density_j);
				sum_value = sum_value - tmp;
			}
		}
		sum_value = sum_value * c;
		Particles[x].m_velocity+=sum_value;
	}
}

void Simulate()
{
	for (unsigned int x = 0; x < Particles.size(); x+=1)
	{
		if(!Particles[x].m_fluids)
			break;

		Particles[x].m_lastPosition=Particles[x].m_position;

		applyGravity(x);
	}

	constraintProjection();
	
	for (unsigned int x = 0; x < Particles.size(); x+=1)
	{
		if(!Particles[x].m_fluids)
			break;

		velocityUpdateFirstOrder(x);
	}

	computeXSPHViscosity();
}

void putSphere(Vector3 pos)
{
    glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(pos.x, pos.y, pos.z);
    glutSolidSphere(rh, 12, 6);
}

void Render()
{
    glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45, 1, 0.1, 200);
	glTranslatef(-15.5, -16, -100);

    for(unsigned int x = 0;x < Particles.size();x += 1)
		if(Particles[x].m_fluids)
			putSphere(Particles[x].m_position);
	
	Simulate();

    glFlush();

	puts("Render Frame Done");
}

void initRigbody()
{
	for (unsigned int x = Particles.size()-1; x >= 0; x-=1)
	{
		if(Particles[x].m_fluids)
			break;
		double delta = poly6WKernel(Vector3(0));
		for (unsigned int y = 0; y < Particles[x].m_neighborhood.size(); y+=1)
		{
			unsigned int neighborIndex = Particles[x].m_neighborhood[y];
			if(Particles[neighborIndex].m_fluids)
				continue;
			delta += poly6WKernel(Particles[x].m_position - Particles[neighborIndex].m_position);
		}
		delta = r0 / delta;
		Particles[x].m_mass=delta;
	}
}

void initParticles()
{
	for(int i=1;i<=fluidWidth;i+=1)
    {
        for(int j=1;j<=fluidWidth;j+=1)
        {
            for(int k=1;k<=fluidWidth;k+=1)
            {
				Particle x;
				x.m_mass=0.8*r0;
				x.m_fluids=true;
				x.m_velocity=Vector3();
				x.m_position=Vector3(10+i,j+20,10+k);
				x.m_neighborhood.clear();
				Particles.push_back(x);
            }
        }
    }
	for(int i=0;i<boxWidth;i+=1)
    {
        for(int j=0;j<boxWidth;j+=1)
        {
            for(int k=0;k<boxWidth;k+=1)
            {
				if(i==0||i==boxWidth-1||j==0||j==boxWidth-1||k==0||k==boxWidth-1)
				{
					Particle x;
					x.m_fluids=false;
					x.m_velocity=Vector3();
					x.m_position=Vector3(i,j,k);
					x.m_neighborhood.clear();
					Particles.push_back(x);
				}
            }
        }
    }
	printf("initParticles Done.\n");
	printf("Particles Count:%d\n",(int)Particles.size());

	updateNeighborhood(false);
	initRigbody();

	printf("Last Wall Mass:%lf\n",Particles[Particles.size()-1].m_mass);
}

void initLight()
{
	GLfloat pos[]{ 30.f,50.f,30.f,0.f };
    GLfloat diff[]{ 0.1f,0.6f,0.7f,0.f };
    GLfloat ambient[]{ 0.2f,0.2f,0.2f,0.f };
    GLfloat specular[]{ 1.f,1.f,1.f,1.f };

    glLightfv(GL_LIGHT0, GL_POSITION, pos);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diff);
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specular);

    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);
}

int main(int argc, char *argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(screenWidth, screenHeight);
	glutCreateWindow("PBF");

    initParticles();
	initLight();

	glutDisplayFunc(&Render);
	glutIdleFunc(&Render);
	glutMainLoop();
	return 0;
}