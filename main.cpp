#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>

#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <GL/glut.h>

using namespace std;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Vector_3 Point_3;
typedef Kernel::Vector_3 Vector_3;

struct CrestPoint {
    Vector_3 point;
    int componentId;
    Vector_3 Ti, To;

};

struct ComposanteConnectee {
    double ridgeness;
    double sphericalness;
    double cyclideness;
};

struct CrestEdge {
    int point1ID,
        point2ID,
        triangleID;
};

struct Spline {
    vector<CrestPoint> points;
};

void affichage(void);

void clavier(unsigned char touche,int x,int y);
void affiche_repere(void);

void mouse(int, int, int, int);
void mouseMotion(int, int);
float t=.5 ;


Spline TCB_Spline(vector<CrestPoint> &P, double tension, double bias, double continuité);
vector<CrestPoint> read();
void drawTCBSpline(Spline courbe, double tension, double continuite, double bias, Vector_3 color);
Vector_3 calculate_To( Vector_3 Pi, Vector_3 Pi_plus_1,Vector_3 Pi_plus_2,double tension, double bias, double continuité);
Vector_3 calculate_Ti(Vector_3 Pi_moins_1, Vector_3 Pi, Vector_3 Pi_plus_1,double tension, double bias, double continuité);

vector<Spline> split_splines(const vector<CrestPoint>& P, int pointID);
double calculate_norme(Vector_3 v);
void similitude(Spline s1, Spline s2);


// variables globales pour OpenGL
bool mouseLeftDown;
bool mouseRightDown;
bool mouseMiddleDown;
bool mouseWheelUp;
bool mouseWheelDown;
float mouseX, mouseY;
float cameraAngleX;
float cameraAngleY;
float cameraDistance=0.;
// constantes pour les materieux
  float no_mat[] = {0.0f, 0.0f, 0.0f, 1.0f};
    float mat_ambient[] = {0.7f, 0.7f, 0.7f, 1.0f};
    float mat_ambient_color[] = {0.8f, 0.8f, 0.2f, 1.0f};
    float mat_diffuse[] = {0.1f, 0.5f, 0.8f, 1.0f};
    float mat_specular[] = {1.0f, 1.0f, 1.0f, 1.0f};
    float no_shininess = 0.0f;
    float low_shininess = 5.0f;
    float high_shininess = 100.0f;
    float mat_emission[] = {0.3f, 0.2f, 0.2f, 0.0f};

Vector_3 red = {1, 0, 0};
Vector_3 blue = {0, 0, 1};
Vector_3 green = {0, 1, 0};



void initOpenGl() 
{ 

//lumiere 

	glClearColor( .5, .5, 0.5, 0.0 );
 
	glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  GLfloat l_pos[] = { 3.,3.5,3.0,1.0 };
  glLightfv(GL_LIGHT0,GL_POSITION,l_pos);

  glLightfv(GL_LIGHT0,GL_DIFFUSE,l_pos);
 glLightfv(GL_LIGHT0,GL_SPECULAR,l_pos);
 glEnable(GL_COLOR_MATERIAL);

  glDepthFunc(GL_LESS);
  glEnable(GL_DEPTH_TEST);
//glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
// glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE|GLUT_RGB);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
         gluPerspective(45.0f,(GLfloat)200/(GLfloat)200,0.1f,10.0f);
	glMatrixMode(GL_MODELVIEW);
      gluLookAt(0.,0.,4., 0.,0.,0., 0.,1.,0.);

}

void similitude(Spline s1, Spline s2)
{
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout <<"EVALUATION SIMILITUDE" << endl;

    auto size1 = s1.points.size(),
         size2 = s2.points.size();


    Vector_3 To = s1.points.at(size1 - 3).To,
             Ti = s2.points.at(1).Ti;
            
    auto normeTi = calculate_norme(Ti);
    auto normeTo = calculate_norme(To); 

    // auto rapportNorme = normeTo / normeTi;

    cout << "Dernier point de la spline 1 : " << s1.points.at(size1-3).point << endl;
    cout << "\tTo = "<< To << " de norme " << normeTo << endl;

    cout << "Premier point de la spline 2 : " << s2.points.at(1).point << endl;
    cout << "\tTi = "<< Ti << " de norme " << normeTi << endl;

    // cout << "\t rapport normeTo / normeTi = " << rapportNorme << endl;

    auto angle = approximate_angle(To, Ti);
    cout << "angle approx = " << round(angle) << endl;
}

double calculate_norme(Vector_3 v) {
    return sqrt(v.squared_length());
}

Vector_3 calculate_Ti(Vector_3 Pi_moins_1, Vector_3 Pi, Vector_3 Pi_plus_1,double tension, double bias, double continuité)
{
    Vector_3 Ti = ((1 - tension) * (1 + bias) * (1 + continuité) * (Pi - Pi_moins_1) +
                       (1 - tension) * (1 - bias) * (1 - continuité) * (Pi_plus_1 - Pi)) / 2.0;
        // cout << "TI = " << Ti << endl;

    return Ti;
}

Vector_3 calculate_To( Vector_3 Pi, Vector_3 Pi_plus_1,Vector_3 Pi_plus_2,double tension, double bias, double continuité)
{
    Vector_3 To = ((1 - tension) * (1 + bias) * (1 - continuité) * (Pi_plus_1 - Pi) +
                       (1 - tension) * (1 - bias) * (1 + continuité) * (Pi_plus_2 - Pi_plus_1)) / 2.0;
    // cout << "TO = " << To << endl;
    return To;
}

vector<Spline> split_splines(const vector<CrestPoint>& P, int pointID)
{
    Spline spline1, spline2;
    vector<Spline> splitSplines;
    if (pointID < 0 || pointID >= P.size()) {
        cout << "in split_splines pointID < 0 || pointID >= P.size()" << endl;
        return splitSplines;
    }
    for (int i = 0; i <= pointID; ++i) {
        // Ajouter les éléments au premier vecteur
        spline1.points.push_back(P[i]);
    }
    for (int i = pointID; i < P.size(); ++i) {
        // Ajouter les éléments au deuxième vecteur
        spline2.points.push_back(P[i]);
    }

    splitSplines.push_back(spline1);
    splitSplines.push_back(spline2);

    return splitSplines;

}


Spline TCB_Spline(vector<CrestPoint> &P, double tension, double bias, double continuité) {
    Spline courbe;
    for (size_t i = 1; i < P.size() - 2; ++i) {
        // Points P_{i-1}, P_i, P_{i+1}
        Vector_3 Pi_moins_1 = P[i - 1].point;
        Vector_3 Pi = P[i].point;
        Vector_3 Pi_plus_1 = P[i + 1].point;

        // Calculer les tangentes Ti et To
        Vector_3 Ti = calculate_Ti(Pi_moins_1, Pi, Pi_plus_1, tension, bias, continuité);
        Vector_3 To = calculate_To(Pi, Pi_plus_1,P[i + 2].point, tension, bias, continuité);

        // Interpoler la courbe entre P_i et P_{i+1}
        for (double t = 0.0; t <= 1.0; t += 0.1) {
            CrestPoint H;
            H.point = (2 * t * t * t - 3 * t * t + 1) * Pi+

                        (t * t * t - 2 * t * t + t) * Ti +

                        (-2 * t * t * t + 3 * t * t) * Pi_plus_1 +

                        (t * t * t - t * t) * To;
            // cout << "Interpolated Point: " << H.point << endl;
            courbe.points.push_back(H);
        }

    }

    auto sizec = courbe.points.size();
    for(int i=1; i<sizec-2; i++)
    {
        // cout << "\t\t\t here !" << endl;
        Vector_3 Ti = calculate_Ti(courbe.points.at(i-1).point, courbe.points.at(i).point, courbe.points.at(i+1).point, tension, bias, continuité);
        Vector_3 To = calculate_To(courbe.points.at(i).point, courbe.points.at(i+1).point, courbe.points.at(i+2).point, tension, bias, continuité);
    
        courbe.points.at(i).Ti = Ti;
        courbe.points.at(i).To = To;
    }

    return courbe;
}

vector<CrestPoint> read() {
    ifstream file("ridges.txt");

    if(!file.is_open()) {
        cerr << "erreur ouverture fichier" << endl;
    }

    //V: nb points de crete
    //E: nb aretes de crete
    //N: nb composantes connectees
    int V, E, N;
    file >> V >> E >> N; //trois premiers entiers du fichiers

    //save crest points in vector crestPoints
    vector<CrestPoint> crestPoints(V);
    for(int i=0; i<V; i++) {
        file >> crestPoints[i].point >> crestPoints[i].componentId;
    }

    //save crest points in vector crestPoints
    vector<ComposanteConnectee> composantesConnectee(N);
    for(int i=0; i<N; i++) {
        file >> composantesConnectee[i].ridgeness;
        file >> composantesConnectee[i].sphericalness;
        file >> composantesConnectee[i].cyclideness;
    }

    //save crest points in vector crestPoints
    vector<CrestEdge> CrestEdges(E);
    for(int i=0; i<E; i++) {
        file >> CrestEdges[i].point1ID;
        file >> CrestEdges[i].point2ID;
        file >> CrestEdges[i].triangleID;
    }

    // cout << "Points de crete :" << endl;
    // for(const auto& cp : crestPoints) {
    //     cout << "Point: " << cp.point << ", ID composantes : " << cp.componentId << endl;
    // } 

    // cout << "Aretes de crete :" << endl;
    // for(const auto& ce : CrestEdges) {
    //     cout << "Edge: " << ce.point1ID << " - " << ce.point2ID << ", ID Triangle : " << ce.triangleID << endl;
    // } 

    // cout << "Composantes connectees :" << endl;
    // for(const auto& cc : composantesConnectee) {
    //     cout << "Composante connectée : Rigidité : " << cc.ridgeness << " Sphéricité : " << cc.sphericalness << " Cyclidicité : " << cc.cyclideness << endl; 
    // } 

    file.close();
    return crestPoints;
}

void drawTCBSpline(Spline courbe, double tension, double continuite, double bias, Vector_3 color) {

    glBegin(GL_LINE_STRIP);

    glColor3f(color.x(), color.y(), color.z());

    for (const CrestPoint& cp : courbe.points) {
        glVertex3f(cp.point.x(), cp.point.y(), cp.point.z());
    }

    glEnd();
}

void afficherTCBSpline()
{
    // vector<CrestPoint> P = {
    //         {Vector_3(-1, -1, 0), 1},
    //         {Vector_3(0, 0, 0.5), 1},
    //         {Vector_3(1, 1, 1), 1},
    //         {Vector_3(1, 2, 1.5), 1}
    //     };

    // vector<CrestPoint> allPoints = {
    //     {Vector_3(-1,-1,-1), 1},
    //     {Vector_3(-0.5,-0.5,-0.5), 1},
    //     {Vector_3(0,0,0), 1},
    //     {Vector_3(0.5,0.5,0.5), 1},
    //     {Vector_3(1,1,1), 1},
    //     {Vector_3(1.5,1.5,1.5), 1}
    // };

    double tension = 0.5;
    double continuity = 0.5;
    double bias = 0.5;

    vector<CrestPoint> allPoints = read();
    if (!allPoints.empty()) {
        CrestPoint dernierPoint = allPoints.back();
        int N = dernierPoint.componentId;
        for(int i=1; i<N+1; i++) {
            vector<CrestPoint> P;
            for(const auto& p:allPoints) {
                if(p.componentId == 7 && P.size() < 20) { //
                    P.push_back(p);
                }
            }

            Spline s = TCB_Spline(P, tension, bias, continuity);
            // drawTCBSpline(s, tension, continuity, bias);
            // cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl << endl << "SPLINE  "  << endl;
            // for(const auto& cp : s.points) {
            //     cout << "Point: " << cp.point << endl << " Ti = " << cp.Ti << " To = " << cp.To << endl;
            // } 
            auto numPoints = s.points.size();
            vector<Spline> splines = split_splines(s.points, numPoints / 2);

            // auto numPoints = P.size();
            // vector<Spline> splines = split_splines(P, numPoints / 2);

            Spline s1 = splines[0];
            Spline s2 = splines[1];
            auto size1 = s1.points.size();
            auto size2 = s2.points.size();

            // for(const auto& cp : s1.points) {
            //     cout << "Point: " << cp.point << endl;
            // } 
            // // cout << "size1 = " << size1 << " and P[size-1] " << s1.points.at(size1-1).point << endl;
            s1 = TCB_Spline(s1.points, tension, bias, continuity);
            s2 = TCB_Spline(s2.points, tension, bias, continuity);

            drawTCBSpline(s1, tension, continuity, bias, red);
            drawTCBSpline(s2, tension, continuity, bias, blue);

            cout << endl << "SPLINE 1 " << endl << "numPoints = " << s1.points.size() << endl;
            
            // for(const auto& cp : s1.points) {
            //     cout << "Point: " << cp.point <<  endl << "Ti = " << cp.Ti << " To = " << cp.To << endl;
            // } 
            
            cout << "SPLINE 2 " << endl << "numPoints = " << s1.points.size()  << endl;
            // for(const auto& cp : s2.points) {
            //     cout << "Point: "  << cp.point << endl << "Ti = " << cp.Ti << " To = " << cp.To << endl;
            // } 

            similitude(s1, s2);
        }

    } else {
        std::cerr << "Le vecteur allPoints est vide." << std::endl;
    }


}

int main(int argc,char **argv) {
  /* initialisation de glut et creation
     de la fenetre */
  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_RGB);
  glutInitWindowPosition(200,200);
  glutInitWindowSize(600,600);
  glutCreateWindow("ifs");

  /* Initialisation d'OpenGL */
  glClearColor(0.0,0.0,0.0,0.0);
  glColor3f(1.0,1.0,1.0);
  glPointSize(1.0);
	
	//ifs = new Ifs();
  /* enregistrement des fonctions de rappel */
  glutDisplayFunc(affichage);
  glutKeyboardFunc(clavier);
  glutMouseFunc(mouse);
  glutMotionFunc(mouseMotion);
  //-------------------------------


  //-------------------------------
    initOpenGl() ;
//-------------------------------

/* Entree dans la boucle principale glut */
  glutMainLoop();
  return 0;
}

void affiche_repere(void)
{
  glBegin(GL_LINES);
  glColor3f(1.0,0.0,0.0);
  glVertex2f(0.,0.);
  glVertex2f(1.,0.);
  glEnd(); 

	 glBegin(GL_LINES);
  glColor3f(0.0,1.0,0.0);
  glVertex2f(0.,0.);
  glVertex2f(0.,1.);
  glEnd(); 
   glBegin(GL_LINES);
  glColor3f(0.0,0.0,1.0);
  glVertex3f(0.,0.,0.);
  glVertex3f(0.,0.,1.);
  glEnd(); 
}

void affichage(void)
{
	glMatrixMode(GL_MODELVIEW);
  /* effacement de l'image avec la couleur de fond */
//	glClear(GL_COLOR_BUFFER_BIT);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//       glClearDepth(10.0f);                         // 0 is near, >0 is far

        glPushMatrix();
	glTranslatef(0,0,cameraDistance);
	glRotatef(cameraAngleX,1.,0.,0.)	;
	glRotatef(cameraAngleY,0.,1.,0.);
	// affiche_repere();
afficherTCBSpline();
        glPopMatrix();
  /* on force l'affichage du resultat */

          glFlush();
  glutSwapBuffers();

}

//------------------------------------------------------


//------------------------------------------------------
void clavier(unsigned char touche,int x,int y)
{

  switch (touche)
    {
    case '+': //
      t+=.1;
       if (t > 1 ) t=1;
      glutPostRedisplay();
      break;
    case '-': //* ajustement du t
       t-=.1;
        if (t < 0 ) t=0;
      glutPostRedisplay();
      break;
    case 'f': //* affichage en mode fil de fer 
      glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
      glutPostRedisplay();
      break;
      case 'p': //* affichage du carre plein 
      glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
      glutPostRedisplay();
      break;
  case 's' : //* Affichage en mode sommets seuls 
      glPolygonMode(GL_FRONT_AND_BACK,GL_POINT);
      glutPostRedisplay();
      break;

    case 'q' : //*la touche 'q' permet de quitter le programme 
      exit(0);
    }
    
}

void mouse(int button, int state, int x, int y)
{
    mouseX = x;
    mouseY = y;

    if(button == GLUT_LEFT_BUTTON)
    {
        if(state == GLUT_DOWN)
        {
            mouseLeftDown = true;
        }
        else if(state == GLUT_UP)
            mouseLeftDown = false;
    }

    else if(button == GLUT_RIGHT_BUTTON)
    {
        if(state == GLUT_DOWN)
        {
            mouseRightDown = true;
        }
        else if(state == GLUT_UP)
            mouseRightDown = false;
    }

    else if(button == GLUT_MIDDLE_BUTTON)
    {
        if(state == GLUT_DOWN)
        {
            mouseMiddleDown = true;
        }
        else if(state == GLUT_UP)
            mouseMiddleDown = false;
    }
    else if (button == 3) // Scroll Up
    {
        if(state == GLUT_DOWN)
        {
            mouseWheelUp = true;
        }
        else if(state == GLUT_UP)
            mouseWheelUp = false;
    }
    else if (button == 4) // Scroll Down
    {
        if(state == GLUT_DOWN)
        {
            mouseWheelDown = true;
        }
        else if(state == GLUT_UP)
            mouseWheelDown = false;
    }
}


void mouseMotion(int x, int y)
{
    if(mouseLeftDown)
    {
        cameraAngleY += (x - mouseX);
        cameraAngleX += (y - mouseY);
        mouseX = x;
        mouseY = y;
    }
    if(mouseRightDown)
    {
        cameraDistance += (y - mouseY) * 0.2f;
        mouseY = y;
    }
    if (mouseWheelUp) // Scroll Up
    {
        cameraDistance -= 1.0f; // Augmentez la valeur si vous voulez un zoom plus rapide
    }
    if (mouseWheelDown) // Scroll Down
    {
        cameraDistance += 1.0f; // Augmentez la valeur si vous voulez un zoom plus rapide
    }

    glutPostRedisplay();
}

    
    
