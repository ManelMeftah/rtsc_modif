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
    int componentId;
    int id;
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

struct Segment {
    CrestPoint P0,
               P1;
    Vector_3 A,
             B,
             C, 
             D;
};


struct Spline {
    vector<CrestPoint> points;
    TCBParameters tcbParams[2];
    vector<Segment> segments;
};


void affichage(void);

void clavier(unsigned char touche,int x,int y);
void affiche_repere(void);

void mouse(int, int, int, int);
void mouseMotion(int, int);
float t=.5 ;

int evaluate_curvature(Segment s1, Segment s2, double seuil) ;

Spline TCB_Spline(vector<CrestPoint> &P, double tension, double continuité, double bias);
vector<CrestPoint> read();
void drawTCBSpline(Spline courbe, double tension, double continuite, double bias, Vector_3 color);
Vector_3 calculate_To( Vector_3 Pi, Vector_3 Pi_plus_1,Vector_3 Pi_moins_1,double tension, double continuité,  double bias);
Vector_3 calculate_Ti(Vector_3 Pi_moins_1, Vector_3 Pi, Vector_3 Pi_plus_1,double tension, double continuité, double bias);

vector<Spline> split_splines(const vector<CrestPoint>& P, int pointID);
double calculate_norme(Vector_3 v);
Vector_3 calculate_vecteur_unitaire(Vector_3 v);

void similitude(Spline s1, Spline s2);
void similitude_iterative(Spline s1, Spline s2);

vector<Vector_3> CalculateSplineSegmentCoefficients( Vector_3 P0, Vector_3 P1,
                                        double delta_k, Vector_3 tangent_o, Vector_3 tangent_i);
vector<pair<int, int>> findEndpoints(const vector<CrestPoint>& crestPoints, const vector<CrestEdge>& crestEdges);
vector<CrestPoint> traceCrestLine(const vector<CrestPoint>& crestPoints, const vector<CrestEdge>& crestEdges, int startingPointID, int componentID);

Vector_3 calculate_courbure(Segment s, double sk);

vector<CrestPoint> read_file();
vector<CrestPoint> read();

double determinant(const Vector_3& u, const Vector_3& v, const Vector_3& w);

double determinant(const Vector_3& u, const Vector_3& v, const Vector_3& w) {
    return CGAL::to_double(CGAL::determinant(u, v, w));
}

Vector_3 round_vector(Vector_3 v);

// variables globales pour OpenGL
bool mouseLeftDown;
bool mouseRightDown;
bool mouseMiddleDown;
bool mouseWheelUp;
bool mouseWheelDown;
float mouseX, mouseY;
float cameraAngleX;
float cameraAngleY;
float cameraDistance=0.1;
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
vector<Segment> splineSegments;

Spline s, s1, s2;

double tension = -1;
double continuity = 0; 
double bias = 0;

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

vector<Vector_3> CalculateSplineSegmentCoefficients(Vector_3 P0,Vector_3 P1,
                                        double delta_k,  Vector_3 tangent_o,  Vector_3 tangent_i) {
    
    vector<Vector_3> coeffs;

    // Vector_3 to = calculate_vecteur_unitaire(tangent_o);
    // Vector_3 ti = calculate_vecteur_unitaire(tangent_i);

    Vector_3 to = tangent_o;
    Vector_3 ti = tangent_i;
    // cout << "\ttangent_o " << to << " delta_k : "<< delta_k << endl;


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
    Vector_3 A = P0;

    // Coefficient B is equal to the product of the time difference Δk and the outgoing tangent Tk^o
    Vector_3 B = to * delta_k;
    // cout << "Vector_3 B = to * delta_k; " << B << endl;

    // Coefficient C is equal to three times the difference between the next control point Pk+1
    // and the current control point Pk, minus the product of the time difference Δk and the sum
    // of twice the outgoing tangent Tk^o and the incoming tangent Tk+1^i
    Vector_3 C = 3 * (P1 - P0) - delta_k * (2 * to + ti);

    // Coefficient D is equal to negative two times the difference between the next control point Pk+1
    // and the current control point Pk, plus the product of the time difference Δk and the sum
    // of the outgoing tangent Tk^o and the incoming tangent Tk+1^i
    Vector_3 D = -2 * (P1 - P0) + delta_k * (to + ti);

    coeffs.push_back(A);
    coeffs.push_back(B);
    coeffs.push_back(C);
    coeffs.push_back(D);
    // cout << "CalculateSplineSegmentCoefficients done" << endl;
    // cout << "\tCalculateSplineSegmentCoefficients." << endl;
    // cout <<"\tA: " << A << " B: "<< B<<" C: " << C << " D: " << D << endl; 

    return coeffs;

}

void display_curvature(Vector_3 curvature) {
    //std::cout << std::fixed << std::setprecision(5); // Réglez la précision décimale selon vos besoins

    std::cout << curvature.x() << " " << curvature.y() << " " << curvature.z() << std::endl;
}


Vector_3 calculate_courbure(Segment s, double sk) {

    // cout << "segment p0 " << s.P0.point << " P1 " << s.P1.point << endl;
    // cout << "A = " << s.A << " B " << s.B << " C " << s.C << " D " << s.D << endl;


    double delta_k = 1.0;

    // Calcul de la première dérivée
//     Vector_3 derivative_first = s.B + 2 * t * s.C + 3 * t2 * s.D;

    Vector_3 first_derivative = 
        1.0 / delta_k * s.B + 2.0 * ((sk - 1) / delta_k) * s.C +
        3.0 * pow(((sk - 1) / delta_k), 2) * s.D;

    // cout << "\tf' " << first_derivative << endl;

    // Calcul de la deuxième dérivée (courbure)
    Vector_3 curvature = 1.0 / pow(delta_k, 2) * (2 * s.C + 6 * ((sk - 1) / delta_k) * s.D);

    return curvature;


    // Vector_3 courbure_moyenne{0, 0, 0};

    // for(double t = 0; t <= 1; t += 0.1) {
    //     double t2 = t * t;

    //     //  C'(t) = B + 2tC + 3t^2D
    //     Vector_3 derivative_first = s.B + 2 * t * s.C + 3 * t2 * s.D;

    //     // C''(t) = 2C + 6tD
    //     Vector_3 derivative_second = 2 * s.C + 6 * t * s.D;
    //     double mixte_product = determinant(derivative_first, derivative_second, Vector_3(0, 0, 0));

    //     // (unused) k(t) = ||C''(t)|| / ||C'(t)||^3 
    //     // k(t) = det(C'(t), C''(t)) / ||C'(t)||^3
    //     double norm_derivative_first = calculate_norme(derivative_first);
    //     double norm_cubed_derivative_first = norm_derivative_first * norm_derivative_first * norm_derivative_first;
    //     // Vector_3 courbure = derivative_second / norm_cubed_derivative_first;
    //     Vector_3 courbure = mixte_product / norm_cubed_derivative_first * derivative_second;
    //     courbure_moyenne += courbure;
    // }
    
    // return courbure_moyenne / 10.0;
}


Vector_3 calculate_vecteur_unitaire(Vector_3 v) 
{
    double norme = calculate_norme(v);
    Vector_3 v_unitaire = v / norme;
    return v_unitaire;
}

Vector_3 round_vector(Vector_3 v) {
    double x = round(v.x());
    double y = round(v.y());
    double z = round(v.z());
    return Vector_3(x, y, z);
}


void similitude(Spline s1, Spline s2) {
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout <<"EVALUATION SIMILITUDE" << endl;

    auto size1 = s1.points.size(),
         size2 = s2.points.size();

    CrestPoint P1 = s1.points.at(size1 - 3),
               P2 = s2.points.at(1);


    Vector_3 To = round_vector(P1.To),
             Ti = round_vector(P2.Ti);
            
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

    cout << "num segs : s1 " << s1.segments.size()<< " s2 " << s2.segments.size() << endl;
    
    int s1_segment_size = s1.segments.size(),
        s2_segment_size = s2.segments.size();

    Segment s1_last_segment = s1.segments.at(s1_segment_size-1),
            s2_first_segment = s2.segments.at(0);
    cout << "last segment of s1 " << s1_last_segment.P0.point << " and " << s1_last_segment.P1.point << endl;
    cout << "first segment of s2 " << s2_first_segment.P0.point << " and " << s2_first_segment.P1.point << endl;

    Vector_3 courbure_s1 = calculate_courbure(s1_last_segment, 1.0);
    Vector_3 courbure_s2 = calculate_courbure(s2_first_segment, 0.0); 

    cout << "\tcoubure s1 = " ;
    display_curvature(courbure_s1);
    cout << "\tcoubure s2 = ";
    display_curvature(courbure_s2);
}

void similitude_iterative(Spline s1, Spline s2) {
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout <<"EVALUATION SIMILITUDE ITERATIVE" << endl;

    auto size1 = s1.points.size();
    auto size2 = s2.points.size();

    // Comparaison entre le dernier point de s1 et le premier point de s2
    CrestPoint P1 = s1.points.at(size1 - 3),
               P2 = s2.points.at(1);

    Vector_3 To = round_vector(P1.To);
    Vector_3 Ti = round_vector(P2.Ti);

    auto normeTi = calculate_norme(Ti);
    auto normeTo = calculate_norme(To);

    auto angle = approximate_angle(To, Ti);
    cout << "Angle approximatif : " << round(angle) << " degrés" << endl;

    // Vector_3 courbure_s1 = calculate_courbure(s1.segments[size1 - 2], 1.0); // Courbure au dernier point de s1
    // Vector_3 courbure_s2 = calculate_courbure(s2.segments[1], 0.0); // Courbure au premier point de s2
    


    int s1_segment_size = s1.segments.size(),
        s2_segment_size = s2.segments.size();


    int i = s1_segment_size-1, j = 0, seuil = 1;
    int score_courbure = 0, score_total = 0;
    bool remaining_points = true;

    do {
        Segment s1_last_segment = s1.segments.at(i),
                s2_first_segment = s2.segments.at(j);

        score_courbure = evaluate_curvature(s1_last_segment, s2_first_segment, seuil);
        cout << "\tscore de similitude courbure : " << score_courbure << endl;
        if(score_courbure >= 5)
            score_total++;
        i--;
        j++;
        seuil++;
        remaining_points = ((i != -1) && (j != s2_segment_size));
    } while(score_courbure >= 5 && remaining_points);
    cout << "score de similitude total : " << score_total << endl;

    // for (double s=0.0; s<=1.0; s = s+0.1)
    // {
    //     Vector_3 courbure_s1 = calculate_courbure(s1_last_segment, 1.0-s);
    //     Vector_3 courbure_s2 = calculate_courbure(s2_first_segment, s); 

    //     // Affichage des courbures
    //     cout << "Point " << 1.0-s << " : ";
    //     cout << "Courbure de s1 : ";
    //     display_curvature(courbure_s1);
        
    //     cout << "Point " << s << " : ";
    //     cout << "Courbure de s2 : ";
    //     display_curvature(courbure_s2);
    //     cout << endl;

    // }
}

int evaluate_curvature(Segment s1, Segment s2, double seuil) {
    // cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    // cout << "EVALUATION DE LA COURBURE ITERATIVE" << endl;

    int score = 0;
    double s = 0.0;
    while (s <= 1.0) {
        // Calcul des courbures au point courant
        Vector_3 courbure_s1 = calculate_courbure(s1, 1.0 - s);
        Vector_3 courbure_s2 = calculate_courbure(s2, s);

        // Affichage des courbures
        cout << "Point " << 1.0 - s << " : " << endl;
        cout << "Courbure de s1 : ";
        display_curvature(courbure_s1);
        cout << "Courbure de s2 : ";
        display_curvature(courbure_s2);
        cout << endl;

        // Comparaison des courbures et mise à jour du score si la différence est inférieure au seuil
        Vector_3 diff = courbure_s1 - courbure_s2;
        double x = round(abs(diff.x()));
        double y = round(abs(diff.y()));
        double z = round(abs(diff.z()));


        Vector_3 difference = Vector_3(x, y, z);   
        cout << " difference : " << difference << endl;

        if (difference.x() <= seuil && difference.y() <= seuil && difference.z() <= seuil) {
            score++;
        } else {
            break; // Sortie de la boucle si la différence dépasse le seuil
        }

        // Passage au point suivant
        s += 0.1;
    }

    return score;
}




double calculate_norme(Vector_3 v) {
    return sqrt(v.squared_length());
}

Vector_3 calculate_Ti(Vector_3 Pi_moins_1, Vector_3 Pi, Vector_3 Pi_plus_1,double tension, double continuité, double bias)
{
        Vector_3 Ti = (1 - tension) * (1 + continuité) * (1 - bias) / 2.0 * (Pi_plus_1 - Pi) +
              (1 - tension) * (1 - continuité) * (1 + bias) / 2.0 * (Pi - Pi_moins_1);

        // Ti = ( Pi_plus_1 - Pi_moins_1 ) / 2.0;
        // cout << "Point: " << Pi << endl;
        // cout << "\tTI = " << Ti << endl;

    return Ti;
}

Vector_3 calculate_To( Vector_3 Pi, Vector_3 Pi_plus_1,Vector_3 Pi_moins_1,double tension, double continuité, double bias)
{
    Vector_3 To = (1 - tension) * (1 - continuité) * (1 - bias) / 2.0 * (Pi_plus_1 - Pi) +
                (1 - tension) * (1 + continuité) * (1 + bias) / 2.0 * (Pi - Pi_moins_1);
    // cout << "\tTO = " << To << endl;
    // To = ( Pi_plus_1 - Pi_moins_1) / 2.0;
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


Spline TCB_Spline(vector<CrestPoint> &P, double tension, double continuité, double bias) {
    Spline courbe;
    for (int i = 0; i < P.size() - 1; i++) {
        // Points P_{i-1}, P_i, P_{i+1}
        Vector_3 Pi_moins_1 = P[i - 1].point;
        Vector_3 Pi = P[i].point;
        Vector_3 Pi_plus_1 = P[i + 1].point;

        // Calculer les tangentes Ti et To
        P[i].Ti = calculate_Ti(Pi_moins_1, Pi, Pi_plus_1, tension, continuité, bias);
        // P[i].To = calculate_To(Pi, Pi_plus_1,P[i + 2].point, tension, continuité ,bias);
        P[i].To = calculate_To(Pi, Pi_plus_1,Pi_moins_1, tension, continuité ,bias);
        cout << "Point: " << Pi << endl;
        cout << "\tTI = " << P[i].Ti << endl;
        cout << "\tTO = " << P[i].To << endl;
        cout << "~~~~~~~~~~~~~~~~~~~" << endl;
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

    Segment segment;
    for (int i = 0; i < P.size() - 1; i++) {
        Vector_3 Pi = P[i].point;
        Vector_3 Pi_plus_1 = P[i + 1].point;



        // cout << "\tpoint n " << i << " and " << i+1 << endl;

        double delta_k = 1.0;
        double sk = 0;
        double sk_plus_1 = 1.0;
        vector<Vector_3> vec_coeffs = CalculateSplineSegmentCoefficients(
            Pi,
            Pi_plus_1,
            delta_k,
            P[i].To,
            P[i + 1].Ti
            
        );

         
        Vector_3 A = vec_coeffs.at(0),
                B = vec_coeffs.at(1),
                C = vec_coeffs.at(2),
                D = vec_coeffs.at(3);


        segment.P0 =  P[i];
        segment.P1 =  P[i+1];
        segment.A = A;
        segment.B = B;
        segment.C = C;
        segment.D = D;

        courbe.segments.push_back(segment);

    

        // cout << "\tA: " << A << " B: " << B << " C: " << C << " D: " << D << endl; 
        // cout << "\ttangent o " << P[i].To << endl;

        for (double s = 0; s <= 1.0; s += 0.1) {
            double t = (s-0)/delta_k; 

            double t2 = t * t;
            double t3 = t2 * t;

            // auto x = A.x() + t * B.x() + t2 * C.x() + t3 * D.x();
            // auto y = A.y() + t * B.y() + t2 * C.y() + t3 * D.y();
            // auto z = A.z() + t * B.z() + t2 * C.z() + t3 * D.z();

            Vector_3 cp =A + t * B + t2 * C + t3 * D;
            if (!std::isnan(cp.x()) && !std::isnan(cp.y()) && !std::isnan(cp.z())){
                CrestPoint H;
                H.point = cp;
                courbe.points.push_back(H);
            }
        }
    }

    auto sizec = courbe.points.size();
    for(int i=1; i<sizec-1; i++)
    {
        // cout << "\t\t\t here !" << endl;
        Vector_3 Ti = calculate_Ti(courbe.points.at(i-1).point, courbe.points.at(i).point, courbe.points.at(i+1).point, tension, continuité, bias );
        Vector_3 To = calculate_To(courbe.points.at(i).point, courbe.points.at(i+1).point, courbe.points.at(i-1).point, tension, continuité, bias);
    
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
// glPointSize(1); // Remplacez "taille_des_points" par la taille souhaitée
    glBegin(GL_LINE_STRIP );

    glColor3f(color.x(), color.y(), color.z());

    for (const CrestPoint& cp : courbe.points) {
        glVertex3f(cp.point.x(), cp.point.y(), cp.point.z());
    }

    glEnd();

}

void drawPoints(vector<CrestPoint> P, Vector_3 color)
{
    
    glColor3f(color.x(), color.y(), color.z());
    for (const CrestPoint& cp : P) {
        glPushMatrix();
        glTranslatef(cp.point.x(), cp.point.y(), cp.point.z());
        glutSolidSphere(0.01, 20, 20); // Dessine une sphère à chaque point avec le rayon spécifié
        glPopMatrix();
    }
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

void initTCBSpline()
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

    
    // vector<CrestPoint> allPoints = {
    //     {Vector_3(-1,-1,-1), 1},
    //     {Vector_3(-0.5,0.0,-0.5), 1},
    //     {Vector_3(0,0,0.5), 1},
    //     {Vector_3(0.5,1,0.5), 1},
    //     {Vector_3(1.5,1,1), 1},
    //     {Vector_3(1.5,2,1.5), 1},
    //     {Vector_3(2,2,1.5), 1}
    // };


    vector<CrestPoint> allPoints = {
        // {Vector_3(-1, 2, 0), 1},
        {Vector_3(-1, 1, 0), 1},
        {Vector_3(0, 0, 0), 1},
        {Vector_3(1, 1, 0), 1}
        // {Vector_3(1, 2, 0), 1}
    };



    // vector<CrestPoint> allPoints = read();
    // vector<CrestPoint> allPoints = read_file();
    if (!allPoints.empty()) {
        CrestPoint dernierPoint = allPoints.back();
        int N = dernierPoint.componentId;
    drawPoints(allPoints, red);

    cout << "N" << N << endl;
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

            //  trouver l'endpoint avec componentID=2.
            for (const auto& endpoint : endpoints) {
                if (endpoint.second == 2) {
                    endpointID = endpoint.first;
                    endpointComponentID = endpoint.second;
                    break; 
                }
            }

            cout << "firstEndpointID " << endpointID << endl;
            cout << "firstEndpointComponentID "  <<  endpointComponentID << endl;

            // vector<CrestPoint> orderedPoints = traceCrestLine(P, crestEdges, endpointID, endpointComponentID);
             vector<CrestPoint> orderedPoints = P;
            s = TCB_Spline(orderedPoints, tension,  continuity, bias);

            cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl << endl << "SPLINE  "  << endl;
            // for(const auto& cp : orderedPoints) {
            //     cout << "Point: " << cp.point << endl ;
            // } 
            // Spline s = TCB_Spline(P, tension, bias, continuity);
    
            // drawTCBSpline(s, tension, bias, continuity, red);
            // cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl << endl << "SPLINE  "  << endl;
            // for(const auto& cp : s.points) {
            //     cout << "Point: " << cp.point << endl << " Ti = " << cp.Ti << " To = " << cp.To << endl;
            // } 
            auto numPoints = s.points.size();
            // drawTCBSpline(s, tension, continuity, bias, blue);

            vector<Spline> splines = split_splines(s.points, numPoints / 2);

            // auto numPoints = P.size();
            // vector<Spline> splines = split_splines(P, numPoints / 2);

            s1 = splines[0];
            s2 = splines[1];
            auto size1 = s1.points.size();
            auto size2 = s2.points.size();

            // for(const auto& cp : s1.points) {
            //     cout << "Point: " << cp.point << endl;
            // } 
            // // cout << "size1 = " << size1 << " and P[size-1] " << s1.points.at(size1-1).point << endl;

            s1 = TCB_Spline(s1.points, tension, bias, continuity);
            s2 = TCB_Spline(s2.points, tension, bias, continuity);


            // cout << endl << "SPLINE 1 " << endl << "numPoints = " << s1.points.size() << endl;
            
            // for(const auto& cp : s1.points) {
            //     cout << "Point: " << cp.point <<  endl << "Ti = " << cp.Ti << " To = " << cp.To << endl;
            // } 
            
            // cout << "SPLINE 2 " << endl << "numPoints = " << s1.points.size()  << endl;
            // for(const auto& cp : s2.points) {
            //     cout << "Point: "  << cp.point << endl << "Ti = " << cp.Ti << " To = " << cp.To << endl;
            // } 

        }


    } else {
        std::cerr << "Le vecteur allPoints est vide." << std::endl;
    }
}

void afficherTCBSpline()
{
    // drawTCBSpline(s, tension, continuity, bias, red);
    
    drawTCBSpline(s1, tension, continuity, bias, red);
    drawTCBSpline(s2, tension, continuity, bias, blue);
    // drawPoints(s1.points, green);
    // drawPoints(s2.points, red);
    drawPoints(s.points, green);

    // Vector_3 courbure = calculate_courbure(s.segments.at(0), 0.0);
    // display_curvature(courbure);
    similitude_iterative(s1, s2);

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
  initTCBSpline();
	
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
	glTranslatef(0,-1,cameraDistance-1);
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

    
    
