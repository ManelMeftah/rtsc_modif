#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stack>

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
    int id;
    int componentId;
    Vector_3 Ti, To;
    //id = line num - 4
};

struct ComposanteConnectee {
    int id;
    double ridgeness;
    double sphericalness;
    double cyclideness;
};

struct CrestEdge {
    int point1ID,
        point2ID,
        triangleID;
};

struct TCBParameters {
    double tension;
    double continuity;
    double bias;
};

struct Spline {
    vector<CrestPoint> points;
    TCBParameters tcbParams[2];
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

vector<Vector_3> CalculateSplineSegmentCoefficients(const Vector_3& P0, const Vector_3& P1,
                                        double delta_k, const Vector_3& tangent_o, const Vector_3& tangent_i);
vector<pair<int, int>> findEndpoints(const vector<CrestPoint>& crestPoints, const vector<CrestEdge>& crestEdges);
vector<CrestPoint> traceCrestLine(const vector<CrestPoint>& crestPoints, const vector<CrestEdge>& crestEdges, int startingPointID, int componentID);

vector<CrestPoint> read_file();
vector<CrestPoint> read();

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

vector<CrestEdge> crestEdges; 

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

vector<Vector_3> CalculateSplineSegmentCoefficients(const Vector_3& P0, const Vector_3& P1,
                                        double delta_k, const Vector_3& tangent_o, const Vector_3& tangent_i) {
    
    vector<Vector_3> coeffs;

    // if (calculate_norme(tangent_o) == 0.0) {
    //     cerr << "Error: Outgoing tangent vector (tangent_o) should not be a zero vector." << endl;
    //     return coeffs;  // Retourner un vecteur vide en cas d'erreur
    // }

    // // Vérification de la non-nullité de la vector tangent_i
    // if (calculate_norme(tangent_i) == 0.0) {
    //     cerr << "Error: Incoming tangent vector (tangent_i) should not be a zero vector." << endl;
    //     return coeffs;  // Retourner un vecteur vide en cas d'erreur
    // }

    
    // cout << "CalculateSplineSegmentCoefficients ..." << endl;
    // Coefficient A is equal to the control point Pk
    Vector_3 A = P1;

    // Coefficient B is equal to the product of the time difference Δk and the outgoing tangent Tk^o
    Vector_3 B = tangent_o * delta_k;

    // Coefficient C is equal to three times the difference between the next control point Pk+1
    // and the current control point Pk, minus the product of the time difference Δk and the sum
    // of twice the outgoing tangent Tk^o and the incoming tangent Tk+1^i
    Vector_3 C = 3 * (P1 - P0) - delta_k * (2 * tangent_o + tangent_i);

    // Coefficient D is equal to negative two times the difference between the next control point Pk+1
    // and the current control point Pk, plus the product of the time difference Δk and the sum
    // of the outgoing tangent Tk^o and the incoming tangent Tk+1^i
    Vector_3 D = -2 * (P1 - P0) + delta_k * (tangent_o + tangent_i);
    coeffs.push_back(A);
    coeffs.push_back(B);
    coeffs.push_back(C);
    coeffs.push_back(D);
    // cout << "CalculateSplineSegmentCoefficients done" << endl;

    return coeffs;

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

    //
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
        P[i].Ti = calculate_Ti(Pi_moins_1, Pi, Pi_plus_1, tension, bias, continuité);
        P[i].To = calculate_To(Pi, Pi_plus_1,P[i + 2].point, tension, bias, continuité);
    }
        // Interpoler la courbe entre P_i et P_{i+1}

/*
        for (double t = 0.0; t <= 1.0; t += 0.1) {
            CrestPoint H;
            H.point = (2 * t * t * t - 3 * t * t + 1) * Pi+

                        (t * t * t - 2 * t * t + t) * Ti +

                        (-2 * t * t * t + 3 * t * t) * Pi_plus_1 +

                        (t * t * t - t * t) * To;
            // cout << "Interpolated Point: " << H.point << endl;
            courbe.points.push_back(H);
        }
*/
    for (size_t i = 1; i < P.size() - 2; ++i) {
        // Points P_{i-1}, P_i, P_{i+1}
        Vector_3 Pi_moins_1 = P[i - 1].point;
        Vector_3 Pi = P[i].point;
        Vector_3 Pi_plus_1 = P[i + 1].point;

        auto sizec = P.size();
        // for (double sk = 1.0; sk <= P.size(); sk += 1.0) {
        // Coefficients pour le segment actuel
            // double sk = static_cast<double>(j);
            // cout << "here : " << sk << endl;

            double delta_k = 1.0;  
            double sk = i,
                   sk_plus_1 = i+1;

            // cout << "\tcourbe.points[i].point " << Pi << " courbe.points[i + 1].point " << Pi_plus_1 << endl;


            vector<Vector_3> coeffs = CalculateSplineSegmentCoefficients(
                Pi,
                Pi_plus_1,
                delta_k, 
                P[i+1].Ti,
                P[i].To
            );

            // courbe.points.at(i).Ti = B;
            // courbe.points.at(i).To = D;
            Vector_3 A = coeffs.at(0),
                     B = coeffs.at(1),
                     C = coeffs.at(2),
                     D = coeffs.at(3);

            // cout << "\tA " << A << " B " << B << " C " << C << " D " << D << endl;

            for (double s = sk; s <= sk_plus_1; s += 0.1) { 
                double t = (s - sk) / delta_k;
                double t2 = t * t;
                double t3 = t2 * t;

                auto x = A.x() + t * B.x() + t2 * C.x() + t3 * D.x();
                auto y = A.y() + t * B.y() + t2 * C.y() + t3 * D.y();
                auto z = A.z() + t * B.z() + t2 * C.z() + t3 * D.z();
                Vector_3 Xk{x, y, z};

                CrestPoint H;
                H.point = Xk;
                if (H.point.x() != 0.0 && H.point.y() != 0.0 && H.point.z() != 0.0) {
                    // cout << "Interpolated Point: " << H.point << endl;
                    courbe.points.push_back(H);
                }
              

            }

        // }
      
/*   */   




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

vector<CrestPoint> read_file() {
    ifstream file("ridges.txt");
        
    if (!file.is_open()) {
        cerr << "erreur ouverture fichier" << endl;
        return {}; // Return an empty vector on error
    }

    // Read header information
    int V, E, N;
    file >> V >> E >> N;


    // Read crest points
    vector<CrestPoint> crestPoints(V);
    for (int i = 0; i < V; i++) {
        file >> crestPoints[i].point >> crestPoints[i].componentId;

        // Assign correct ID based on line number:
        crestPoints[i].id = i; 
    }

    
    // Read connected components
    vector<ComposanteConnectee> composantesConnectee(N);
    for (int i = 0; i < N; i++) {
        file >> composantesConnectee[i].ridgeness;
        file >> composantesConnectee[i].sphericalness;
        file >> composantesConnectee[i].cyclideness;
        composantesConnectee[i].id = 0;
    }

    // Read crest edges
    // vector<CrestEdge> crestEdges(E);
    crestEdges.resize(E);
    for (int i = 0; i < E; i++) {
        file >> crestEdges[i].point1ID;
        file >> crestEdges[i].point2ID;
        file >> crestEdges[i].triangleID;
    }


    file.close();

    // // Handle connecting edges (-1 triangleID)
    // for (auto& edge : crestEdges) {
    //     if (edge.triangleID == -1) {
    //         // This is a connecting edge, adjust point IDs if needed
    //         // (Logic for adjusting IDs based on component connectivity)
    //     }
    // }

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

vector<pair<int, int>> findEndpoints(const vector<CrestPoint>& crestPoints, const vector<CrestEdge>& crestEdges) {
    vector<pair<int, int>> endpoints;
    unordered_map<int, int> neighborCounts; // Map point ID to neighbor count

    for (const CrestEdge& edge : crestEdges) {
        neighborCounts[edge.point1ID]++;
        neighborCounts[edge.point2ID]++;
    }

    for (int i = 0; i < crestPoints.size(); i++) {
        if (neighborCounts[crestPoints[i].id] == 1) {
            endpoints.push_back({crestPoints[i].id, crestPoints[i].componentId});
        }
    }

    return endpoints;
}

vector<CrestPoint> traceCrestLine(const vector<CrestPoint>& crestPoints, const vector<CrestEdge>& crestEdges, int startingPointID, int componentID) {
    vector<CrestPoint> orderedPoints;
    stack<int> pointStack;
    unordered_set<int> visitedPoints;

    pointStack.push(startingPointID);

    while (!pointStack.empty()) {
        int currentPointID = pointStack.top();
        pointStack.pop();

        orderedPoints.push_back(crestPoints[currentPointID]);
        visitedPoints.insert(currentPointID);

        for (const CrestEdge& edge : crestEdges) {
            if (edge.point1ID == currentPointID && visitedPoints.count(edge.point2ID) == 0 && crestPoints[edge.point2ID].componentId == componentID) {
                pointStack.push(edge.point2ID);
            } else if (edge.point2ID == currentPointID && visitedPoints.count(edge.point1ID) == 0 && crestPoints[edge.point1ID].componentId == componentID) {
                pointStack.push(edge.point1ID);
            }
        }
    }

    return orderedPoints;
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

    // vector<CrestPoint> allPoints = read();
    vector<CrestPoint> allPoints = read_file();
    if (!allPoints.empty()) {
        CrestPoint dernierPoint = allPoints.back();
        int N = dernierPoint.componentId;


    for(int i=1; i<N+1; i++) {
            vector<CrestPoint> P;
            for(const auto& p:allPoints) {
                // if(p.componentId == 2) { // && P.size() < 20
                    P.push_back(p);
                // }
            }

            vector<pair<int, int>> endpoints = findEndpoints(P, crestEdges);
            // int firstEndpointID = endpoints[3].first;
            // int firstEndpointComponentID = endpoints[3].second;

            int endpointID = -1;  // Initialisez à une valeur qui ne peut pas être un ID valide.
            int endpointComponentID = -1;

            // Parcourez le vecteur endpoints pour trouver l'endpoint avec componentID=2.
            for (const auto& endpoint : endpoints) {
                if (endpoint.second == 2) {
                    endpointID = endpoint.first;
                    endpointComponentID = endpoint.second;
                    break;  // Sortez de la boucle dès que vous avez trouvé l'endpoint recherché.
                }
            }

            cout << "firstEndpointID " << endpointID << endl;
            cout << "firstEndpointComponentID "  <<  endpointComponentID << endl;

            vector<CrestPoint> orderedPoints = traceCrestLine(P, crestEdges, endpointID, endpointComponentID);
            Spline s = TCB_Spline(orderedPoints, tension, bias, continuity);

            cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl << endl << "SPLINE  "  << endl;
            for(const auto& cp : orderedPoints) {
                cout << "Point: " << cp.point << endl ;
            } 
    /**
     *  
            Spline s = TCB_Spline(P, tension, bias, continuity);
    */
            // drawTCBSpline(s, tension, continuity, bias);
            // cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl << endl << "SPLINE  "  << endl;
            // for(const auto& cp : s.points) {
            //     cout << "Point: " << cp.point << endl << " Ti = " << cp.Ti << " To = " << cp.To << endl;
            // } 
            auto numPoints = s.points.size();
            // drawTCBSpline(s, tension, continuity, bias, blue);

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

            // cout << endl << "SPLINE 1 " << endl << "numPoints = " << s1.points.size() << endl;
            
            // for(const auto& cp : s1.points) {
            //     cout << "Point: " << cp.point <<  endl << "Ti = " << cp.Ti << " To = " << cp.To << endl;
            // } 
            
            // cout << "SPLINE 2 " << endl << "numPoints = " << s1.points.size()  << endl;
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

    
    
