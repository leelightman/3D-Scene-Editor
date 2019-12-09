// This example is heavily based on the tutorial at https://open.gl
// C++ include
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <list>

// OpenGL Helpers to reduce the clutter
#include "Helpers.h"

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
// GLFW is necessary to handle the OpenGL context
#include <GLFW/glfw3.h>
#else
// GLFW is necessary to handle the OpenGL context
#include <GLFW/glfw3.h>
#endif


// Linear Algebra Library
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>

// Shortcut to avoid  and std::everywhere, DO NOT USE IN .h
using namespace std;
using namespace Eigen;
struct Vertex{
    float x;
    float y;
    float z;
};

struct Face{
    int v_0;
    int v_1;
    int v_2;
};

struct Model{
    Array<Vertex,502,1> V;
    Array<Face,1000,1> F;
    Matrix<Vertex,Dynamic,Dynamic> F_V;
    int vertex_num;
    int face_num;
    Array<Vector3d,1000,1> face_normal;
    Array<Vector3d, 502,1> vertex_normal;
    int model_enum;
    Vector3d center;
    
    Matrix4f trans_mat;
    Matrix4f scale_mat;
    Matrix4f rotate_mat;
    Matrix4f model_mat;
};

// VertexBufferObject wrapper for cube
// VBO_F_pos is the VBO for vertex corresponding to each face
VertexBufferObject VBO,VBO_N,VBO_F_V,VBO_N2;
VertexBufferObject VBO_bunny,VBO_bunny_N,VBO_bunny_F_V,VBO_bunny_N2;
VertexBufferObject VBO_bumpy,VBO_bumpy_N,VBO_bumpy_F_V,VBO_bumpy_N2;
ElementBufferObject EBO,EBO_bunny,EBO_bumpy;

// Contains the vertex positions for cube, 3 rows and 36 columns, each column is a vertex
MatrixXf V_cube(3,8),V_cube_f_v(3,36),V_bunny(3,502),V_bunny_f_v(3,3000),V_bumpy(3,502),V_bumpy_f_v(3,3000);
MatrixXi elements_cube(3,12),elements_bunny(3,1000),elements_bumpy(3,1000);
MatrixXf normal_cube(3,8),normal2_cube(3,36),normal_bunny(3,502),normal2_bunny(3,3000),normal_bumpy(3,502),normal2_bumpy(3,3000);
Vector2i is_selected(3,0);
// Contains the view transformation
Matrix4f view(4,4),projection(4,4);
// Geometry
int cube_ct=0,bunny_ct=0,bumpy_ct=0;
enum rendering_mode{wireframe, flat,phong};
enum model_type{cube_tri,bumpy,bunny};
enum camera_mode_enum {ortho,pers,trackball};
int mode = 1,camera_mode = 0; //init value not 0/1/2

// Object Model
Model cube_model,bunny_model,bumpy_model;
Model cube_models[50],bunny_models[50], bumpy_models[50];
Model empty_model;
// const definition
const double EPSILON = 0.0000001;
const double PI = 3.1415926;
// Camera
Vector3d cameraPos(0.0f, 0.0f, 3.0f),cameraTarget(0.0f, 0.0f, 0.0f),up (0.0f, 1.0f, 0.0f); 
GLfloat camX,camY,camZ;
float xchange=0,ychange=0;
// helper functions
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}
Vector3d tri_normal(Vector3d v1, Vector3d v2, Vector3d v3){
    Vector3d triangle_normal;
    Vector3d edge1 = v2-v1;
    Vector3d edge2 = v3-v1;
    triangle_normal = edge1.cross(edge2).normalized();
    return triangle_normal;
}

void face_normal(Model &instance,MatrixXf &V_N2){

    V_N2.resize(3,instance.face_num*3);
    for(int i = 0; i<instance.face_num;i++){
        int face_index_v1 = instance.F[i].v_0;
        int face_index_v2 = instance.F[i].v_1;
        int face_index_v3 = instance.F[i].v_2;
        Vertex v1_vertex = instance.V[face_index_v1];
        Vertex v2_vertex = instance.V[face_index_v2];
        Vertex v3_vertex = instance.V[face_index_v3];
        Vector3d v1(v1_vertex.x,v1_vertex.y,v1_vertex.z);
        Vector3d v2(v2_vertex.x,v2_vertex.y,v2_vertex.z);
        Vector3d v3(v3_vertex.x,v3_vertex.y,v3_vertex.z);
        instance.face_normal[i] = tri_normal(v1,v2,v3);
        // Add face normal to normal of each point in V_N2
        V_N2.col(3*i)<< instance.face_normal[i][0],instance.face_normal[i][1],instance.face_normal[i][2];
        V_N2.col(3*i+1)=V_N2.col(3*i);
        V_N2.col(3*i+2)=V_N2.col(3*i);
    }
}

list<int> list_of_faces_with_vertex(int index,Model model){
    list<int> ans;
    for(int i = 0; i<model.face_num;i++){
        if(model.F[i].v_0 == index ||model.F[i].v_1 == index||model.F[i].v_2 == index)
            ans.push_back(i);
    }

    //cout <<"face contain index "<<index<<" is ";
    // for (list<int>::iterator j=ans.begin(); j!=ans.end(); j++){
    //     cout << *j << " ";
    // }
    return ans;
}

Vector3d calc_vertex_normal(Model model,list<int> ans){
    Vector3d sum_vec3(0,0,0);
    int counter=0;
    for (list<int>::iterator i=ans.begin(); i!=ans.end(); i++){
        //counter++;
        // step 2: get the normal of each face
        //cout<< *i<<" ";
        sum_vec3+=model.face_normal[*i];
        
    }
    //cout<<"sum_vec3 "<<sum_vec3<<endl;
    Vector3d vec_ans(sum_vec3[0],sum_vec3[1],sum_vec3[2]);
    //cout<<"vec_ans "<<vec_ans.normalized()<<endl;
    return vec_ans.normalized();
}
void vertex_normal(Model &model,MatrixXf &V_N){
    V_N.resize(3,model.vertex_num);
    for(int i = 0; i<model.vertex_num;i++){
        list<int> ans = list_of_faces_with_vertex(i,model);
        Vector3d ver_nor=calc_vertex_normal(model,ans);
        V_N.col(i)<<ver_nor[0],ver_nor[1],ver_nor[2];
    }
    //cout<<"V_N.cols()"<<V_N.col(501)<<endl;
}

Vector3d model_barycenter(Model model){
    Vector3d bc;
    float xsum=0,ysum=0,zsum=0;
    for (int i=0;i<model.vertex_num;i++){
        xsum+=model.V[i].x;
        ysum+=model.V[i].y;
        zsum+=model.V[i].z;
    }
    bc[0] = xsum/model.vertex_num;
    bc[1] = ysum/model.vertex_num;
    bc[2] = zsum/model.vertex_num;
    return bc;
}

void translate_model(Model &model, double xlen,double ylen,double zlen){
    Matrix4f translate_mat;
    translate_mat<<1,    0,    0, xlen,
                   0,    1,    0, ylen,
                   0,    0,    1, zlen,
                   0,    0,    0,    1;
    model.model_mat*=translate_mat;
}

void scale_model(Model &model, double scale){
    Matrix4f scale_mat;
    scale_mat<<scale,   0,      0,  0,
               0,   scale,      0,  0,
               0,       0,  scale,  0,
               0,       0,      0,  1; 
    model.model_mat*=scale_mat;
}

void rotate_model(Model &model,int degree){
    float angle = degree*M_PI/180;
    //translate_model(model,-model.center[0],-model.center[1],-model.center[2]);
    Matrix4f rotation_mat;
    rotation_mat<< 1,          0,          0,0,
                   0, cos(angle), sin(angle),0,
                   0,-sin(angle), cos(angle),0,
                   0,          0,          0,1;
    model.model_mat=model.model_mat*rotation_mat;
    //translate_model(model,model.center[0],model.center[1],model.center[2]);
}

void delete_model(Model &model){
    model = empty_model;
}


void calculate_world_position (Vector3d cursor_pos, Vector3f &world_pos, double width, double height){
    // Convert screen position to world coordinates
    
    Vector4f p_screen(cursor_pos[0],height-1-cursor_pos[1],-1.0,1);
    Vector4f p_canonical((p_screen[0]/width)*2-1,(p_screen[1]/height)*2-1,0,1);
    Vector4f p_world = view.inverse()*p_canonical;
    world_pos[0] = p_world[0];
    world_pos[1] = p_world[1];
    world_pos[2] = p_world[2];
}
// Calculate cursor's position in camera space
void calculate_cam_position (Vector3d cursor_pos, Vector3f &cur_cam_pos, double width, double height){
    // Convert screen position to world coordinates
    
    Vector4f p_screen(cursor_pos[0],height-1-cursor_pos[1],-1.0,1);
    Vector4f p_canonical((p_screen[0]/width)*2-1,(p_screen[1]/height)*2-1,-1,1);
    //Vector4f p_world = view.inverse()*p_canonical;
    cur_cam_pos[0] = p_canonical[0];
    cur_cam_pos[1] = p_canonical[1];
    cur_cam_pos[2] = p_canonical[2];
}

bool ray_tri_intersect(Vector3f ray_origin, Vector3f ray_direction, VectorXf v0,VectorXf v1, VectorXf v2, double &t){
    Vector3f edge1, edge2, h, s, q;
    double a,f,u,v;
    edge1 = v1 - v0;
    edge2 = v2 - v0;
    h = ray_direction.cross(edge2);
    a = edge1.dot(h);
    if (a > -EPSILON && a < EPSILON)
        return false;  
    f = 1.0/a;
    s = ray_origin - v0;
    u = f * s.dot(h);
    if (u < 0.0 || u > 1.0)
        return false;
    q = s.cross(edge1);
    v = f * ray_direction.dot(q);
    if (v < 0.0 || u + v > 1.0)
        return false;
    t = f * edge2.dot(q);
    if (t > EPSILON){
        //ray_intersecction = ray_origin + ray_direction * t;
        return true;
    } else
        return false;

}
Vector2i on_object(Vector3f cur_cam_pos){
    // calculate if world_pos is on the object. if mouse is on certain object,
    Vector3f ray_dir;
    if (camera_mode==ortho) {
        ray_dir={0,0,-1};
    } else if (camera_mode==pers){
        ray_dir = cur_cam_pos.normalized();
    }
    
    VectorXf v0,v1,v2;
    Vector2i ans(3,0);
    double t=0.0;
    double t_temp = 10;
    //cout<<"p4"<<endl;
    for(int i=0;i<cube_ct;i++){
        for (int j=0;j<cube_model.face_num; j++){
            v0.conservativeResize(4);
            v0[0]=V_cube_f_v(0,3*j);
            v0[1]=V_cube_f_v(1,3*j);
            v0[2]=V_cube_f_v(2,3*j);
            v0[3]=1;
            v0 = view*cube_models[i].model_mat*v0;
            v0.conservativeResize(3);
            //cout<<"p2"<<endl;
            v1.conservativeResize(4);
            v1[0] = V_cube_f_v.col(3*j+1)[0];
            v1[1] = V_cube_f_v.col(3*j+1)[1];
            v1[2] = V_cube_f_v.col(3*j+1)[2];
            v1[3]=1;
            v1= view*cube_models[i].model_mat*v1;
            v1.conservativeResize(3);
            
            v2.conservativeResize(4);
            v2[0] = V_cube_f_v.col(3*j+2)[0];
            v2[1] = V_cube_f_v.col(3*j+2)[1];
            v2[2] = V_cube_f_v.col(3*j+2)[2];
            v2[3]=1;
            v2= view*cube_models[i].model_mat*v2;
            v2.conservativeResize(3);
            bool answer = ray_tri_intersect(cur_cam_pos, ray_dir, v0, v1, v2, t);
            //cout<<"point1: "<<answer<<endl;
            if (answer and t<t_temp){
                ans[0] = cube_tri;
                ans[1] = i;
                t_temp = t;
                //cout<<"cube ans : "<<ans<<endl;
                //return ans;
            }
        }
    }
    for(int i=0;i<bunny_ct;i++){
        for (int j=0;j<bunny_model.face_num; j++){
            v0.conservativeResize(4);
            v0[0]=V_bunny_f_v(0,3*j);
            v0[1]=V_bunny_f_v(1,3*j);
            v0[2]=V_bunny_f_v(2,3*j);
            v0[3]=1;
            v0 = view*bunny_models[i].model_mat*v0;
            v0.conservativeResize(3);
            //cout<<"p2"<<endl;
            v1.conservativeResize(4);
            v1[0] = V_bunny_f_v.col(3*j+1)[0];
            v1[1] = V_bunny_f_v.col(3*j+1)[1];
            v1[2] = V_bunny_f_v.col(3*j+1)[2];
            v1[3]=1;
            v1= view*bunny_models[i].model_mat*v1;
            v1.conservativeResize(3);
            
            v2.conservativeResize(4);
            v2[0] = V_bunny_f_v.col(3*j+2)[0];
            v2[1] = V_bunny_f_v.col(3*j+2)[1];
            v2[2] = V_bunny_f_v.col(3*j+2)[2];
            v2[3]=1;
            v2= view*bunny_models[i].model_mat*v2;
            v2.conservativeResize(3);
            bool answer = ray_tri_intersect(cur_cam_pos, ray_dir, v0, v1, v2, t);
            //cout<<"point1: "<<answer<<endl;
            if (answer and t<t_temp){
                ans[0] = bunny;
                ans[1] = i;
                t_temp = t;
                //cout<<ans<<endl;
                //return ans;
            }
        }
    }
    for(int i=0;i<bumpy_ct;i++){
        for (int j=0;j<bumpy_model.face_num; j++){
            v0.conservativeResize(4);
            v0[0]=V_bumpy_f_v(0,3*j);
            v0[1]=V_bumpy_f_v(1,3*j);
            v0[2]=V_bumpy_f_v(2,3*j);
            v0[3]=1;
            v0 = view*bumpy_models[i].model_mat*v0;
            v0.conservativeResize(3);
            //cout<<"p2"<<endl;
            v1.conservativeResize(4);
            v1[0] = V_bumpy_f_v.col(3*j+1)[0];
            v1[1] = V_bumpy_f_v.col(3*j+1)[1];
            v1[2] = V_bumpy_f_v.col(3*j+1)[2];
            v1[3]=1;
            v1= view*bumpy_models[i].model_mat*v1;
            v1.conservativeResize(3);
            
            v2.conservativeResize(4);
            v2[0] = V_bumpy_f_v.col(3*j+2)[0];
            v2[1] = V_bumpy_f_v.col(3*j+2)[1];
            v2[2] = V_bumpy_f_v.col(3*j+2)[2];
            v2[3]=1;
            v2= view*bumpy_models[i].model_mat*v2;
            v2.conservativeResize(3);
            bool answer = ray_tri_intersect(cur_cam_pos, ray_dir, v0, v1, v2, t);
            //cout<<"point1: "<<answer<<endl;
            if (answer and t< t_temp){
                ans[0] = bumpy;
                ans[1] = i;
                t_temp = t;
                //cout<<ans<<endl;
                //return ans;
            }   
        }
    }
    return ans;
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    // Get the position of the mouse in the window

    double xpos, ypos, xworld, yworld;
    Vector3d cursor_pos;
    Vector3f world_pos;
    Vector3f cur_cam_pos;
    // Get the size of the window
    int width, height;
    glfwGetWindowSize(window, &width, &height);
    glfwGetCursorPos(window, &cursor_pos[0], &cursor_pos[1]);
    calculate_cam_position (cursor_pos, cur_cam_pos, width, height);
    //calculate_world_position (cursor_pos, world_pos, width, height);
    // Update the position of the first vertex if the left button is pressed
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS){
        //cout<<"p3"<<endl;
        //is_selected = on_object(world_pos);
        is_selected = on_object(cur_cam_pos);
        //cout<<"is_selected[0]: "<<is_selected[0]<<endl;
        if (is_selected[0]==cube_tri){
            cout<<"Cube "<< is_selected[1] <<" is selected"<<endl;
            mode = flat;
        } else if (is_selected[0]==bunny){
            cout<<"Bunny "<< is_selected[1] <<" is selected"<<endl;
            mode = flat;
        } else if (is_selected[0]==bumpy){
            cout<<"Bumpy "<< is_selected[1] <<" is selected"<<endl;
            mode = flat;
        }
    }
}

/*
 * Read OFF file format
 * 
 * @param filename
 * */
void read_off_file(string filename, Model &model_instance){
    string readLine;
    ifstream inFile;
    inFile.open(filename);
    if(!inFile){
        cerr << "Unable to open"<<filename<<endl;
        exit(1);   // call system to stop
    } else{
        getline(inFile, readLine);
        getline(inFile, readLine);
        
        // Find out number of vertices
        int space1 = readLine.find(" ", 0);
        int vertex_number = atoi(readLine.substr(0,space1+1).c_str());
        
        // Find out number of faces
        int space2 = readLine.find(" ", space1+1);
        int face_number = atoi(readLine.substr(space1,space2+1).c_str());
        int space3 = 0;
        int space4 = 0;

        // Read vertices into array V
        Array<Vertex, 502,1> V;
        for (int n=0; n<vertex_number; n++){
            getline(inFile, readLine);
            space1 = readLine.find(" ", 0); 
            space2 = readLine.find(" ", space1+1);
            space3 = readLine.find("\n", space2+1);

            V[n].x = atof(readLine.substr(0,space1+1).c_str());
            V[n].y = atof(readLine.substr(space1,space2+1).c_str());
            V[n].z = atof(readLine.substr(space2).c_str());
        }

        // read faces into array F
        Array<Face, 1000,1> F;
        Matrix<Vertex,3,1000> F_V;
        for (int n=0; n<face_number; n++){
            getline(inFile, readLine);
            space1 = readLine.find(" ", 0);
            space2 = readLine.find(" ", space1+1);
            space3 = readLine.find(" ", space2+1);
            space4 = readLine.find(" ", space3+1);

            F[n].v_0 = atoi(readLine.substr(space1,space2+1 ).c_str());            
            F[n].v_1 = atoi(readLine.substr(space2,space3+1 ).c_str());            
            F[n].v_2 = atoi(readLine.substr(space3).c_str());
    
            F_V.col(n)[0]=V[F[n].v_0];
            F_V.col(n)[1]=V[F[n].v_1];
            F_V.col(n)[2]=V[F[n].v_2];
        }
        model_instance.V = V;
        model_instance.F = F;
        model_instance.vertex_num = vertex_number;
        model_instance.face_num = face_number;
        model_instance.F_V = F_V;
        model_instance.center = model_barycenter(model_instance);
        model_instance.model_mat << 1,0, 0, 0,
                                    0,1, 0, 0,
                                    0,0, 1, 0,
                                    0,0, 0, 1;
    }
    inFile.close();
}
void load_VBO(string filename,MatrixXf &V_obj,MatrixXf &V_obj_f_v, MatrixXi &elements_obj){
    string readLine;
    ifstream inFile;
    inFile.open(filename);
    if(!inFile){
        cerr << "Unable to open"<<filename<<endl;
        exit(1);   // call system to stop
    } else{
        getline(inFile, readLine);
        getline(inFile, readLine);
        
        // Find out number of vertices
        int space1 = readLine.find(" ", 0);
        int vertex_number = atoi(readLine.substr(0,space1+1).c_str());
        
        // Find out number of faces
        int space2 = readLine.find(" ", space1+1);
        int face_number = atoi(readLine.substr(space1,space2+1).c_str());
        int space3 = 0;
        int space4 = 0;
        
        V_obj.resize(3,vertex_number);
        elements_obj.resize(3,face_number);
        V_obj_f_v.resize(3,face_number*3);

        // Read vertices into array V
        
        for (int n=0; n<vertex_number; n++){
            getline(inFile, readLine);
            space1 = readLine.find(" ", 0); 
            space2 = readLine.find(" ", space1+1);
            space3 = readLine.find("\n", space2+1);
            float x = atof(readLine.substr(0,space1+1).c_str());
            float y = atof(readLine.substr(space1,space2+1).c_str());
            float z = atof(readLine.substr(space2).c_str());
            V_obj.col(n)<<x,y,z; 
        }

        for (int n=0; n<face_number; n++){
            getline(inFile, readLine);
            space1 = readLine.find(" ", 0);
            space2 = readLine.find(" ", space1+1);
            space3 = readLine.find(" ", space2+1);
            space4 = readLine.find(" ", space3+1);

            int v_0 = atoi(readLine.substr(space1,space2+1 ).c_str());            
            int v_1 = atoi(readLine.substr(space2,space3+1 ).c_str());            
            int v_2 = atoi(readLine.substr(space3).c_str());
            elements_obj.col(n)<<v_0,v_1,v_2;

            V_obj_f_v.col(3*n)<<V_obj.col(v_0);
            V_obj_f_v.col(3*n+1)<<V_obj.col(v_1);
            V_obj_f_v.col(3*n+2)<<V_obj.col(v_2);
        }
    }
    inFile.close();
}

void load_model(){
    cube_model.model_enum = cube_tri;
    read_off_file("../data/cube_tri.off",cube_model);
    translate_model(cube_model,-0.5,-0.5,-0.5);
    
    face_normal(cube_model,normal2_cube);
    vertex_normal(cube_model,normal_cube);
    load_VBO("../data/cube_tri.off",V_cube,V_cube_f_v,elements_cube);

    bumpy_model.model_enum = bumpy;
    read_off_file("../data/bumpy_cube.off",bumpy_model);
    scale_model(bumpy_model, 0.1);

    face_normal(bumpy_model,normal2_bumpy);
    vertex_normal(bumpy_model,normal_bumpy);
    load_VBO("../data/bumpy_cube.off",V_bumpy,V_bumpy_f_v,elements_bumpy);
    

    bunny_model.model_enum = bunny;
    read_off_file("../data/bunny.off",bunny_model);
    scale_model(bunny_model, 4);
    translate_model(bunny_model, -bunny_model.center[0],-bunny_model.center[1],-bunny_model.center[2]);
    
    face_normal(bunny_model,normal2_bunny);
    vertex_normal(bunny_model,normal_bunny);
    load_VBO("../data/bunny.off",V_bunny,V_bunny_f_v,elements_bunny);
}

void add_object(int object_enum){
    switch (object_enum)
    {
        case cube_tri:
            cube_ct +=1;
            cube_models[cube_ct-1] = cube_model;
            break;
        case bumpy:      
            bumpy_ct+=1;
            bumpy_models[bumpy_ct-1] = bumpy_model;
            break;
        case bunny:
            bunny_ct+=1;
            bunny_models[bunny_ct-1] = bunny_model;
            break;
        default:
            break;
    }
}
//void reset_cam_pos(){

//}
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    // Update the position of the first vertex if the keys 1,2, or 3 are pressed
    switch (key)
    {
        case  GLFW_KEY_1:
            if (action == GLFW_PRESS){
                add_object(cube_tri);
                break;
            }
        case  GLFW_KEY_2:
            if (action == GLFW_PRESS){
                add_object(bumpy);
                break;
            }
        case  GLFW_KEY_3:
            if (action == GLFW_PRESS){
                add_object(bunny);
                break;
            }
        case  GLFW_KEY_T:
            if (action == GLFW_PRESS){
                cout<<"is_selected[0] in T: "<<is_selected[0]<<endl;
                float x=0;
                if(mods ==0){
                    x = 0.1;
                } else{
                    x = -0.1;
                }
                switch(is_selected[0]){
                    case cube_tri:
                        translate_model(cube_models[is_selected[1]], x,0,0);
                        break;
                    case bunny:
                        translate_model(bunny_models[is_selected[1]], x*0.1,0,0);
                        break;
                    case bumpy:
                        translate_model(bumpy_models[is_selected[1]], 2*x,0,0);
                        break;
                    default:
                        break;
                }
                cout<<"Translate selected object"<<endl;
            }
            break;
        case  GLFW_KEY_D:
            if (action == GLFW_PRESS){
                cout<<"is_selected[0] in T: "<<is_selected[0]<<endl;
                float x=0;
                switch(is_selected[0]){
                    case cube_tri:
                        delete_model(cube_models[is_selected[1]]);
                        break;
                    case bunny:
                        delete_model(bunny_models[is_selected[1]]);
                        break;
                    case bumpy:
                        delete_model(bumpy_models[is_selected[1]]);
                        break;
                    default:
                        break;
                }
                cout<<"Translate selected object"<<endl;
            }
            break;
        case  GLFW_KEY_S:
            if (action == GLFW_PRESS){
                //scale_model(model, 2);
                //off_to_mat(model,V_bunny);
                float scale = 1;
                if(mods == 0){
                    scale = 0.9;
                } else{
                    scale =1.1;
                }
                Vector3d bc;
                switch(is_selected[0]){
                    case cube_tri:
                        bc = model_barycenter(cube_models[is_selected[1]]);
                        translate_model(cube_models[is_selected[1]], bc[0],bc[1],bc[2]);
                        scale_model(cube_models[is_selected[1]], scale);
                        translate_model(cube_models[is_selected[1]], -bc[0],-bc[1],-bc[2]);
                        break;
                    case bunny:
                        bc = model_barycenter(bunny_models[is_selected[1]]);
                        translate_model(bunny_models[is_selected[1]], bc[0],bc[1],bc[2]);
                        scale_model(bunny_models[is_selected[1]], scale);
                        translate_model(bunny_models[is_selected[1]], -bc[0],-bc[1],-bc[2]);
                        break;
                    case bumpy:
                        scale_model(bumpy_models[is_selected[1]], scale);
                        break;
                    default:
                        break;
                }
                cout<<"Scale selected object"<<endl;
            }
            break;
        case  GLFW_KEY_R:
            if (action == GLFW_PRESS){
                int degree = 15;
                Vector3d bc;
                switch(is_selected[0]){
                    case cube_tri:
                        bc = model_barycenter(cube_models[is_selected[1]]);
                        translate_model(cube_models[is_selected[1]], bc[0],bc[1],bc[2]);
                        rotate_model(cube_models[is_selected[1]], degree);
                        translate_model(cube_models[is_selected[1]], -bc[0],-bc[1],-bc[2]);
                        break;
                    case bunny:
                        bc = model_barycenter(bunny_models[is_selected[1]]);
                        translate_model(bunny_models[is_selected[1]], bc[0],bc[1],bc[2]);
                        rotate_model(bunny_models[is_selected[1]], degree);
                        translate_model(bunny_models[is_selected[1]], -bc[0],-bc[1],-bc[2]);
                        break;
                    case bumpy:
                        rotate_model(bumpy_models[is_selected[1]], degree);
                        break;
                    default:
                        break;
                }
                cout<<"Rotate selected object"<<endl;
            }
            break;
        case  GLFW_KEY_W:
            if (action == GLFW_PRESS){
                mode = wireframe;
                cout<<"Wirefram Mode"<<endl;
            }
            break;
        case  GLFW_KEY_F:
            if (action == GLFW_PRESS){
                mode = flat;
                cout<<"Flat Shading Mode"<<endl;
            }
            break;
        case  GLFW_KEY_P:
            if (action == GLFW_PRESS){
                mode = phong;
                cout<<"Phong Shading Mode"<<endl;
            }
            break;
        case  GLFW_KEY_9:
            if (action == GLFW_PRESS){
                camera_mode = pers;
                cout<<"Perspective camera"<<endl;
            }
            break;
        case  GLFW_KEY_8:
            if (action == GLFW_PRESS){
                camera_mode = trackball;
                //reset_cam_pos();
                cout<<"Trackball camera"<<endl;
            }
            break;
        case  GLFW_KEY_0:
            if (action == GLFW_PRESS){
                camera_mode = ortho;
                cout<<"Orthographic camera"<<endl;
            }
            break;
        case  GLFW_KEY_UP:
            if(action == GLFW_PRESS){
                if(camera_mode == pers || camera_mode == ortho){
                    cameraPos[1] +=0.1;
                } else {
                    ychange += 0.2f;
                }
            }
            break;
        case  GLFW_KEY_LEFT:
            if(action == GLFW_PRESS){
                if(camera_mode == pers || camera_mode == ortho){
                    cameraPos[0] -=0.1;
                } else {
                    xchange -= 0.2f;
                }
            }
            break;
        case  GLFW_KEY_DOWN:
            if(action == GLFW_PRESS) {
                if(camera_mode == pers || camera_mode == ortho){
                    cameraPos[1] -=0.1;
                } else {
                    ychange -= 0.2f;
                }
            }
            break;
        case  GLFW_KEY_RIGHT:
            if(action == GLFW_PRESS){
                if(camera_mode == pers || camera_mode == ortho){
                    cameraPos[0] +=0.1;
                } else {
                    xchange += 0.2f;
                }
            }
            break;
        case GLFW_KEY_ESCAPE:
            if (action == GLFW_PRESS)
            glfwSetWindowShouldClose(window, GL_TRUE);
        default:
            break;
    }
}

int main(void)
{
    GLFWwindow* window;
    double r=1.0,l=-1.0,t=1.0,b=-1.0,n=1,f=1000;

    // Initialize the library
    if (!glfwInit())
        return -1;

    // Activate supersampling
    glfwWindowHint(GLFW_SAMPLES, 8);

    // Require the OpenGL context to support OpenGL 3.2 at the least
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);

    // On apple we have to load a core profile with forward compatibility
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
    // Create a windowed mode window and its OpenGL context
    window = glfwCreateWindow(640, 480, "Hello World", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    // Make the window's context current
    glfwMakeContextCurrent(window);

    #ifndef __APPLE__
      glewExperimental = true;
      GLenum err = glewInit();
      if(GLEW_OK != err)
      {
        /* Problem: glewInit failed, something is seriously wrong. */
       fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
      }
      glGetError(); // pull and savely ignonre unhandled errors like GL_INVALID_ENUM
      fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));
    #endif

    int major, minor, rev;
    major = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MAJOR);
    minor = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MINOR);
    rev = glfwGetWindowAttrib(window, GLFW_CONTEXT_REVISION);
    printf("OpenGL version recieved: %d.%d.%d\n", major, minor, rev);
    printf("Supported OpenGL is %s\n", (const char*)glGetString(GL_VERSION));
    printf("Supported GLSL is %s\n", (const char*)glGetString(GL_SHADING_LANGUAGE_VERSION));
  

    load_model();
    
    // Initialize the VAO
    // A Vertex Array Object (or VAO) is an object that describes how the vertex
    // attributes are stored in a Vertex Buffer Object (or VBO). This means that
    // the VAO is not the actual object storing the vertex data,
    // but the descriptor of the vertex data.
    VertexArrayObject VAO;
    VAO.init();
    VAO.bind();
    VBO.init(); // Initialize the VBO with the vertices data
    VBO_F_V.init();
    EBO.init();
    VBO_N.init();
    VBO_N2.init();  
    VBO_N.update(normal_cube);
    VBO_N2.update(normal2_cube);  
    VBO.update(V_cube); // A VBO is a data container that lives in the GPU memory
    VBO_F_V.update(V_cube_f_v);
    EBO.update(elements_cube);
 
    VertexArrayObject VAO_bunny;
    VAO_bunny.init();
    VAO_bunny.bind();
    VBO_bunny.init();
    VBO_bunny_F_V.init();
    EBO_bunny.init();
    VBO_bunny_N.init();
    VBO_bunny_N2.init();  
    VBO_bunny_N.update(normal_bunny);
    VBO_bunny_N2.update(normal2_bunny); 
    VBO_bunny.update(V_bunny);
    VBO_bunny_F_V.update(V_bunny_f_v);
    EBO_bunny.update(elements_bunny);
 
    VertexArrayObject VAO_bumpy;
    VAO_bumpy.init();
    VAO_bumpy.bind();
    VBO_bumpy.init();
    VBO_bumpy_F_V.init();
    EBO_bumpy.init();
    VBO_bumpy_N.init();
    VBO_bumpy_N2.init();  
    VBO_bumpy_N.update(normal_bumpy);
    VBO_bumpy_N2.update(normal2_bumpy); 
    VBO_bumpy.update(V_bumpy); 
    VBO_bumpy_F_V.update(V_bumpy_f_v);   
    EBO_bumpy.update(elements_bumpy);

    /*
    Initialize the OpenGL Program
    A program controls the OpenGL pipeline and it must contains
    at least a vertex shader and a fragment shader to be valid */
    Program program;

    const GLchar* vertex_shader =
            "#version 150 core\n"
                    "in vec3 position;"
                    "uniform mat4 model;"
                    "uniform mat4 view;"
                    "uniform mat4 proj;"
                    "uniform mat4 camera;"
                    "in vec3 face_normal;"
                    "out vec3 FragPos;"
                    "out vec3 Normal;"
                    "void main()"
                    "{"
                    "    gl_Position = proj*view*model* vec4(position, 1.0);"
                    "    Normal = vec3(model *vec4(face_normal,1.0));"
                    "    FragPos = vec3(model * vec4(position, 1.0f));"
                    "}";
    const GLchar* fragment_shader =
            "#version 150 core\n"
                    "out vec4 outColor;"
                    "in vec3 Normal;"
                    "uniform vec3 lightColor;"
                    "uniform vec3 viewPos;"
                    "uniform vec3 lightPos; "
                    "in vec3 FragPos;"
                    "uniform vec3 color;"
                    "void main()"
                    "{"
                    "    float ambientStrength = 0.2f;"
                    "    vec3 ambient = ambientStrength * lightColor;"
                    "    vec3 norm = normalize(Normal);"
                    "    vec3 lightDir = normalize(lightPos - FragPos);"
                    "    float diff = max(dot(norm, lightDir), 0.0);"
                    "    vec3 diffuse = diff * lightColor;"
                    "    float specularStrength = 0.5f;"
                    "    vec3 viewDir = normalize(viewPos - FragPos);"
                    "    vec3 reflectDir = reflect(-lightDir, norm);"
                    "    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 100);"
                    "    vec3 specular = spec * lightColor;"
                    "    vec3 result = (ambient + diffuse+specular) * color;"
                    "    outColor = vec4(result, 1.0f);"
                    "}";
    
    // Compile the two shaders and upload the binary to the GPU
    // Note that we have to explicitly specify that the output "slot" called outColor
    // is the one that we want in the fragment buffer (and thus on screen)
    // Called helper function
    program.init(vertex_shader,fragment_shader,"outColor");
    program.bind();
    
    glUniform3f(program.uniform("lightPos"), 0.0f,0.0f,1.0f);  
    glUniform3f(program.uniform("lightColor"), 1.0f,1.0f,1.0f);  
    // The vertex shader wants the position of the vertices as an input.
    // The following line connects the VBO we defined above with the position "slot"
    // in the vertex shader

    // Register the keyboard callback
    glfwSetKeyCallback(window, key_callback);
    
    // Register the mouse callback
    //glfwSetCursorPosCallback(window, cursor_position_callback);

    // Register the mouse callback
    glfwSetMouseButtonCallback(window, mouse_button_callback);

    // Update viewport
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glEnable(GL_DEPTH_TEST);

    //Camera
    Vector3d cameraDirection,cameraRight,cameraUp;
    Matrix4f mat1(4,4),mat2(4,4);

/*
    Vector3d cameraPos(0.0f, 0.0f, 3.0f);
    Matrix4f proj_multi(4,4);
    Vector3d cameraTarget(0.0f, 0.0f, 0.0f);
    Vector3d cameraDirection= (cameraTarget - cameraPos).normalized();
    cameraDirection = -cameraDirection;
    Vector3d up (0.0f, 1.0f, 0.0f); 
    Vector3d cameraRight = up.cross(cameraDirection).normalized();
    Vector3d cameraUp= cameraDirection.cross(cameraRight);
    Matrix4f mat1(4,4),mat2(4,4);
    mat1<<
        cameraRight[0],cameraRight[1],cameraRight[2],0,
        cameraUp[0],cameraUp[1],cameraUp[2],0,
        cameraDirection[0],cameraDirection[1],cameraDirection[2],0,
        0,0,0,1;
    mat2<<
        1,0,0,-cameraPos[0],
        0,1,0,-cameraPos[1],
        0,0,1,-cameraPos[2],
        0,0,0,1;
    view <<  mat1*mat2;
*/
    /*
     * Loop until the user closes the window
     * Main loop 
     */

    while (!glfwWindowShouldClose(window))
    {
        // Bind your program, only once
        program.bind();
        
        // Get size of the window
        int width, height;
        glfwGetWindowSize(window, &width, &height);
        float aspect_ratio = float(height)/float(width); // corresponds to the necessary width scaling
        Matrix4f proj_multi(4,4);
        proj_multi <<
            aspect_ratio,0, 0, 0,
            0,           1, 0, 0,
            0,           0, 1, 0,
            0,           0, 0, 1;
        switch(camera_mode){
            case ortho:
                projection <<
                    2/(r-l),0, 0, -(r+l)/(r-l),
                    0,2/(t-b), 0, -(t+b)/(t-b),
                    0,0,-2/(f-n), -(f+n)/(f-n),
                    0,0, 0, 1;
                projection*=proj_multi;
                break;
            case pers:
                projection <<
                    2*n/(r-l), 0, (r+l)/(r-l),0,
                    0,2*n/(t-b), (t+b)/(t-b),0,
                    0,0,-(f+n)/(f-n),-2*f*n/(f-n), 
                    0,0, -1, 0;
                projection*=proj_multi;
                break;
        }
        
        switch(camera_mode){
            case ortho:
                cameraDirection= (cameraTarget - cameraPos).normalized();
                cameraDirection = -cameraDirection;
                cameraRight = up.cross(cameraDirection).normalized();
                cameraUp= cameraDirection.cross(cameraRight);
                mat1<<
                    cameraRight[0],cameraRight[1],cameraRight[2],0,
                    cameraUp[0],cameraUp[1],cameraUp[2],0,
                    cameraDirection[0],cameraDirection[1],cameraDirection[2],0,
                    0,0,0,1;
                mat2<<
                    1,0,0,-cameraPos[0],
                    0,1,0,-cameraPos[1],
                    0,0,1,-cameraPos[2],
                    0,0,0,1;
                view <<  mat1*mat2;
                break;
            case pers:
                cameraDirection= (cameraTarget - cameraPos).normalized();
                cameraDirection = -cameraDirection;
                cameraRight = up.cross(cameraDirection).normalized();
                cameraUp= cameraDirection.cross(cameraRight);
                mat1<<
                    cameraRight[0],cameraRight[1],cameraRight[2],0,
                    cameraUp[0],cameraUp[1],cameraUp[2],0,
                    cameraDirection[0],cameraDirection[1],cameraDirection[2],0,
                    0,0,0,1;
                mat2<<
                    1,0,0,-cameraPos[0],
                    0,1,0,-cameraPos[1],
                    0,0,1,-cameraPos[2],
                    0,0,0,1;
                view <<  mat1*mat2;
                break;
            case trackball:
                GLfloat radius = 10.0f;
                camX = sin(xchange) * radius;
                camY = sin(ychange) * radius;
                camZ = cos(xchange)*cos(ychange) * radius;
                //camX = sin(PI*xchange)* radius;
                //camY = sin(PI*ychange)*radius;
                //camZ = cos(PI*xchange) * radius;
                cameraPos = {camX,camY,camZ};
                cameraDirection= (cameraTarget - cameraPos).normalized();
                cameraDirection = -cameraDirection;
                cameraRight = up.cross(cameraDirection).normalized();
                cameraUp= cameraDirection.cross(cameraRight);
                mat1<<
                    cameraRight[0],cameraRight[1],cameraRight[2],0,
                    cameraUp[0],cameraUp[1],cameraUp[2],0,
                    cameraDirection[0],cameraDirection[1],cameraDirection[2],0,
                    0,0,0,1;
                mat2<<
                    1,0,0,-cameraPos[0],
                    0,1,0,-cameraPos[1],
                    0,0,1,-cameraPos[2],
                    0,0,0,1;
                view <<  mat1*mat2;
                break;
        }
        

        glUniformMatrix4fv(program.uniform("view"), 1, GL_FALSE, view.data()); 
        glUniformMatrix4fv(program.uniform("proj"), 1, GL_FALSE, projection.data()); 
        glUniform3f(program.uniform("viewPos"), cameraPos[0],cameraPos[1],cameraPos[2]);
        // Clear the framebuffer. Set background color
        glClearColor(.5f, .5f, .5f, 1.0f);
        //glClear(GL_COLOR_BUFFER_BIT);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Draw cube
        VAO.bind();
        for(int i=0;i<cube_ct;i++){           
            if(is_selected[0]==cube_tri && is_selected[1] == i){
                glUniform3f(program.uniform("color"), 0.0f,0.0f,1.0f); 
                glUniformMatrix4fv(program.uniform("model"), 1, GL_FALSE, cube_models[i].model_mat.data());
                if(mode == wireframe){
                    program.bindVertexAttribArray("position",VBO_F_V);
                    program.bindVertexAttribArray("face_normal", VBO_N2);
                    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                    glDrawArrays(GL_TRIANGLES,0,V_cube_f_v.cols());
                    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                }else if (mode == flat){
                    program.bindVertexAttribArray("position",VBO_F_V);
                    program.bindVertexAttribArray("face_normal", VBO_N2);
                    glUniform3f(program.uniform("color"), 0.0f,0.0f,0.0f); 
                    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                    glDrawArrays(GL_TRIANGLES,0,V_cube_f_v.cols());
                    glUniform3f(program.uniform("color"), 0.0f,0.0f,1.0f);
                    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                    glDrawArrays(GL_TRIANGLES,0,V_cube_f_v.cols());
                }else if ( mode == phong){
                    program.bindVertexAttribArray("position",VBO);
                    program.bindVertexAttribArray("face_normal", VBO_N);
                    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                    glDrawElements(GL_TRIANGLES,elements_cube.size(),GL_UNSIGNED_INT,0);
                }
            }else{
                program.bindVertexAttribArray("position",VBO_F_V);
                program.bindVertexAttribArray("face_normal", VBO_N2);
                glUniformMatrix4fv(program.uniform("model"), 1, GL_FALSE, cube_models[i].model_mat.data());
                glUniform3f(program.uniform("color"), 0.0f,0.5f,1.0f); 
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                glDrawArrays(GL_TRIANGLES,0,V_cube_f_v.cols());
            }
        }

        // Draw bunny
        VAO_bunny.bind();
        for(int i=0;i<bunny_ct;i++){           
            if(is_selected[0]==bunny && is_selected[1] == i){
                glUniform3f(program.uniform("color"), 0.0f,1.0f,0.0f); 
                glUniformMatrix4fv(program.uniform("model"), 1, GL_FALSE, bunny_models[i].model_mat.data());
                if(mode == wireframe){
                    program.bindVertexAttribArray("position",VBO_bunny_F_V);
                    program.bindVertexAttribArray("face_normal", VBO_bunny_N2);
                    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                    glDrawArrays(GL_TRIANGLES,0,V_bunny_f_v.cols());
                    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                }else if (mode == flat){
                    program.bindVertexAttribArray("position",VBO_bunny_F_V);
                    program.bindVertexAttribArray("face_normal", VBO_bunny_N2);
                    glUniform3f(program.uniform("color"), 0.0f,0.0f,0.0f); 
                    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                    glDrawArrays(GL_TRIANGLES,0,V_bunny_f_v.cols());
                    glUniform3f(program.uniform("color"), 0.0f,1.0f,0.0f);
                    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                    glDrawArrays(GL_TRIANGLES,0,V_bunny_f_v.cols());
                }else if ( mode == phong){
                    program.bindVertexAttribArray("position",VBO_bunny);
                    program.bindVertexAttribArray("face_normal", VBO_bunny_N);
                    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                    glDrawElements(GL_TRIANGLES,elements_bunny.size(),GL_UNSIGNED_INT,0);
                }
            }else{
                program.bindVertexAttribArray("position",VBO_bunny_F_V);
                program.bindVertexAttribArray("face_normal", VBO_bunny_N2);
                glUniformMatrix4fv(program.uniform("model"), 1, GL_FALSE, bunny_models[i].model_mat.data());
                glUniform3f(program.uniform("color"), .5f,.5f,0.0f); 
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                glDrawArrays(GL_TRIANGLES,0,V_bunny_f_v.cols());
            }
        }

        // Draw bumpy
        VAO_bumpy.bind();
        for(int i=0;i<bumpy_ct;i++){           
            if(is_selected[0]==bumpy && is_selected[1] == i){
                glUniform3f(program.uniform("color"), 1.0f,0.0f,0.0f); 
                glUniformMatrix4fv(program.uniform("model"), 1, GL_FALSE, bumpy_models[i].model_mat.data());
                if(mode == wireframe){
                    program.bindVertexAttribArray("position",VBO_bumpy_F_V);
                    program.bindVertexAttribArray("face_normal", VBO_bumpy_N2);
                    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                    glDrawArrays(GL_TRIANGLES,0,V_bumpy_f_v.cols());
                    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                }else if (mode == flat){
                    program.bindVertexAttribArray("position",VBO_bumpy_F_V);
                    program.bindVertexAttribArray("face_normal", VBO_bumpy_N2);
                    glUniform3f(program.uniform("color"), 0.0f,0.0f,0.0f); 
                    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                    glDrawArrays(GL_TRIANGLES,0,V_bumpy_f_v.cols());
                    glUniform3f(program.uniform("color"), 1.0f,0.0f,0.0f);
                    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                    glDrawArrays(GL_TRIANGLES,0,V_bumpy_f_v.cols());
                }else if ( mode == phong){
                    program.bindVertexAttribArray("position",VBO_bumpy);
                    program.bindVertexAttribArray("face_normal", VBO_bumpy_N);
                    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                    glDrawElements(GL_TRIANGLES,elements_bumpy.size(),GL_UNSIGNED_INT,0);
                }
            }else{
                program.bindVertexAttribArray("position",VBO_bumpy_F_V);
                program.bindVertexAttribArray("face_normal", VBO_bumpy_N2);
                glUniformMatrix4fv(program.uniform("model"), 1, GL_FALSE, bumpy_models[i].model_mat.data());
                glUniform3f(program.uniform("color"), 1.0f,0.5f,0.5f); 
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                glDrawArrays(GL_TRIANGLES,0,V_bunny_f_v.cols());
            }
        }
        glfwSwapBuffers(window);
        // Poll for and process events
        glfwPollEvents();
    }

    // Deallocate opengl memory
    program.free();
    VAO.free();
    VAO_bumpy.free();
    VAO_bunny.free();
    VBO.free();
    VBO_bunny.free();
    VBO_bumpy.free();
    VBO_F_V.free();
    VBO_bunny_F_V.free();
    VBO_bumpy_F_V.free();
    VBO_N.free();
    VBO_bunny_N.free();
    VBO_bumpy_N.free();
    VBO_N2.free();
    VBO_bunny_N2.free();
    VBO_bumpy_N2.free();
    EBO.free();
    EBO_bunny.free();
    EBO_bumpy.free();
    //Deallocate glfw internals
    glfwTerminate();
    return 0;
}
