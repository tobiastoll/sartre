//==============================================================================
//  nucleusAnimation.cpp
//
//  Copyright (C) 2010-2019 Tobias Toll and Thomas Ullrich 
//
//  This file is part of Sartre. 
//
//  This program is free software: you can redistribute it and/or modify 
//  it under the terms of the GNU General Public License as published by 
//  the Free Software Foundation.   
//  This program is distributed in the hope that it will be useful, 
//  but without any warranty; without even the implied warranty of 
//  merchantability or fitness for a particular purpose. See the 
//  GNU General Public License for more details. 
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//  Author: Tobias Toll
//  Last update: 
//  $Date: 2019-03-08 20:13:19 +0100 (Fri, 08 Mar 2019) $
//  $Author: ullrich $
//==============================================================================
#include <iostream> 
#include <cmath> 
#include <vector>
#include "TVector3.h"       
#include <algorithm>
#include <unistd.h>
#include <cstdlib>
#include "Nucleus.h"
#include "TableGeneratorNucleus.h"

#include <glut.h>

#define PR(x) cout << #x << " = " << (x) << endl; 


using namespace std; 


// actual vector representing the camera's direction
float lx=0.0f,lz=-1.0f;
// XZ position of the camera
float x=0.0f, y=40.0f, z=5.0f;
// the key states. These variables will be zero
// when no key is being presses
float deltalr=0.0f, deltaud=0.0f;
bool newNucleus=0;
float deltazoom=0.0f;

TableGeneratorNucleus* myNucleus;// = new TableGeneratorNucleus;

void changeSize(int w, int h) {
    
	// Prevent a divide by zero, when window is too short
	// (you cant make a window of zero width).
	if (h == 0)
		h = 1;
	float ratio =  w * 1.0 / h;
    
	// Use the Projection Matrix
	glMatrixMode(GL_PROJECTION);
    
	// Reset Matrix
	glLoadIdentity();
    
	// Set the viewport to be the entire window
	glViewport(0, 0, w, h);
    
	// Set the correct perspective.
	gluPerspective(45.0f, ratio, 0.1f, 100.0f);
    
	// Get Back to the Modelview
	glMatrixMode(GL_MODELVIEW);
}

void drawNucleon() {
    
    //  glColor3f(1.0f, .0f, .0f);
    glTranslatef(0.0f ,0.75f, 0.0f);
    glutSolidSphere(0.7f,20,20);
    
}
 
void computeDir(){
    
    double sud=sin(deltaud), cud=cos(deltaud),
    slr=sin(deltalr), clr=cos(deltalr);
    double zoom=1.+deltazoom;
    //  theta=lr, phi=ud, psi=0;
    x= (clr*x + sud*slr*y + cud*slr*z)*zoom;
    y= (            cud*y +  (-sud)*z)*zoom;
    z=(-slr*x + sud*clr*y + cud*clr*z)*zoom;
}

// Position the light at the origin.
const GLfloat light_pos[] = { 0.0, 50.0, 0.0, 1.0 };
// Attenuation factors for light.
GLfloat const_att = 1.0;

void renderScene(void) {
    
	if (deltalr or deltaud or deltazoom)
        computeDir();
	if( newNucleus )
        while(!myNucleus->generate()){};
    
	// Clear Color and Depth Buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
	// Reset transformations
	glLoadIdentity();
	// Set the camera
	gluLookAt(x, y, z, // camera position
              0., 0.,  0., //camera focus
              0.0f, 1.,  0.0f); //Vector pointing up
	
	for(unsigned int i=0; i<myNucleus->configuration.size(); i++){
        TVector3 pos=myNucleus->configuration.at(i).position();
        glPushMatrix();
        glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, const_att);
        glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
        glTranslatef(pos.X(),pos.Y(),pos.Z());
        drawNucleon();
        glPopMatrix();
	}
	glutSwapBuffers();
}
void idle(void){
    glutPostRedisplay();
}


void pressKey(int key, int xx, int yy) {
    xx=yy;
    switch (key) {
        case GLUT_KEY_LEFT : deltalr = -0.01; break;
        case GLUT_KEY_RIGHT : deltalr = 0.01; break;
        case GLUT_KEY_UP : deltaud = 0.01; break;
        case GLUT_KEY_DOWN : deltaud = -0.01; break;
        case 27 :
        case 81 :
        case 113: exit(0); break;
        case 110: 
        case 78: newNucleus=true; break;
        case 122: 
        case 90: deltazoom = 0.01; break;
        case 65:
        case 97: deltazoom = -0.01; break;
    }
}

void releaseKey(int key, int x, int y) {
    x=y;
    switch (key) {
        case GLUT_KEY_LEFT :
        case GLUT_KEY_RIGHT : deltalr = 0.0; break;
        case GLUT_KEY_UP :
        case GLUT_KEY_DOWN : deltaud = 0; break;
        case 110: 
        case 78: newNucleus=false; break;
        case 122: 
        case 90:
        case 65:
        case 97: deltazoom = 0; break;
    }
}

void exitFunction(){
    delete myNucleus;
}

void usage(const char* prog)
{
    cout << "Usage: " << prog << "[-w] A" << endl;
    cout << "       Valid A are: 208, 197, 110, 63, 40, 27, 16, 1" << endl;
    cout << "       -w    White instead of black background" << endl;
}

int main(int argc, char **argv) {
    
    atexit(exitFunction);
    
    if (argc < 2) {
         usage(argv[0]);
         return 2;
    }
    
    bool hasWhiteBackground = false;
    
    int ch;
    while ((ch = getopt(argc, argv, "w")) != -1) {
        switch (ch) {
            case 'w':
                hasWhiteBackground = true;
                break;
            default:
                usage(argv[0]);
                return 2;
                break;
        }
    }

    if (optind == argc) {
        usage(argv[0]);
        return 2;
    }
    
    string str(argv[optind]);
    for (unsigned int i=0; i<str.length(); i++) {
        if(!isdigit(str[i])) {
            cout << "Error, argument must be a digit." << endl;
            return 2;
        }
    }
    
    unsigned int A = atoi(argv[optind]);
    myNucleus = new TableGeneratorNucleus(A);
    while(!myNucleus->generate()){};
    
    cout<<"List of Commands:"<<endl;
    cout<<"\trotate view: left/right up/down arrow keys"<<endl;
    cout<<"\tzoom in/out: a/z"<<endl;
    cout<<"\tregenerate nucleus: n"<<endl;
    cout<<"\tquit: q/esc"<<endl;
    
    // init GLUT and create window
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100,100);
    glutInitWindowSize(320,320);
    glutCreateWindow("Sartre Nucleus");
    glClearColor(.0f, .0f, .0f, 0.0f); 

    
    //define colours:
    
    //light
    const GLfloat white[]= { 1.0, 1.0, 1.0, 1.0 };
    const GLfloat magenta[]= { 1.0, 0.0, 1.0, 1.0 };
    const GLfloat gold[] = { 1.0, 1.0, 0.0, 1.0 };
    const GLfloat red[] = { 1.0, 0.0, 0.0, 1.0 };
    const GLfloat green[] = { 0.0, 1.0, 0.0, 1.0 };
    const GLfloat blue[] = { 0.0, 0.0, 1.0, 1.0 };
    
    //  const GLfloat color[];
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    switch(A){
        case 208:
        case 16:
            glLightfv(GL_LIGHT0, GL_DIFFUSE, blue);
            glLightfv(GL_LIGHT0, GL_SPECULAR, blue);
            break;
        case 197:
            glLightfv(GL_LIGHT0, GL_DIFFUSE, gold);
            glLightfv(GL_LIGHT0, GL_SPECULAR, gold);
            break;
        case 63:
            glLightfv(GL_LIGHT0, GL_DIFFUSE, red);
            glLightfv(GL_LIGHT0, GL_SPECULAR, red);
            break;
        case 27:
        case 40:
            glLightfv(GL_LIGHT0, GL_DIFFUSE, white);
            glLightfv(GL_LIGHT0, GL_SPECULAR, white);
            break;
        case 110:
            glLightfv(GL_LIGHT0, GL_DIFFUSE, magenta);
            glLightfv(GL_LIGHT0, GL_SPECULAR, magenta);
            break;
        case 1:
            glLightfv(GL_LIGHT0, GL_DIFFUSE, green);
            glLightfv(GL_LIGHT0, GL_SPECULAR, green);
            break;
            
    }
    
    if (hasWhiteBackground)
        glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
    
    // register callbacks
    glutDisplayFunc(renderScene);
    glutReshapeFunc(changeSize);
    glutIdleFunc(idle);
    
    glutSpecialFunc(pressKey);
    
    // here are the new entries
    glutIgnoreKeyRepeat(1);
    glutSpecialUpFunc(releaseKey);
    
    // OpenGL init
    glEnable(GL_DEPTH_TEST);
    
    // enter GLUT event processing cycle
    glutMainLoop();
    
    return 0;
}
