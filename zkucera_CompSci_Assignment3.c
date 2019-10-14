#include <windows.h>
#include <math.h>
#include "camera.c"
#include "fillPoly.c"
const char g_szClassName[] = "myWindowClass";

WNDCLASSEX wc;
MSG Msg;
HWND hwnd;
HDC hdc;
PAINTSTRUCT ps;

#define Lx 3.0 //
#define Ly 5.0 //  Light source coords
#define Lz 3.0 //
#define Ls 1.0 //Intensity of the light source

#define Ia 0.5//Ambient light source

#define Pd 0.50 //Coeff for diffuse light
#define Pa 0.05 //Coeff for ambient light
#define Ps 0.45 //Coeff for specular light

struct polygon {
    dmatrix_t world_points[4];
    dmatrix_t camera_points[4];
    int RED;
    int GREEN;
    int BLUE;
    dmatrix_t normal;
    dmatrix_t centroid;
    dmatrix_t s;
    dmatrix_t r;
    dmatrix_t v;
    double Id;
    double Is;
    float distanceFromCamera;
}p;


//Module Name: swap
//Author: https://www.geeksforgeeks.org/quick-sort/
//Date: March 12th, 2019
//Purpose: Used to swap two elements
//Parameters: a,b: the two elements to be swapped
void swap(struct polygon* a, struct polygon* b){ 
    struct polygon t = *a; 
    *a = *b; 
    *b = t; 
} 

//Module Name: partition
//Author: https://www.geeksforgeeks.org/quick-sort/
//Date: March 12th, 2019
//Purpose: Performs the "pivot" of quick sort
//Parameters: arr[]: the array to be sorted, low: the lowest point of the sort, high: the max size of the array
int partition (struct polygon arr[], int low, int high) { 
    float pivot = arr[high].distanceFromCamera;    // pivot 
    int i = (low - 1);  // Index of smaller element 
  
    for (int j = low; j <= high- 1; j++) 
    { 
        // If current element is smaller than or 
        // equal to pivot 
        if (arr[j].distanceFromCamera >= pivot) 
        { 
            i++;    // increment index of smaller element 
            swap(&arr[i], &arr[j]); 
        } 
    } 
    swap(&arr[i + 1], &arr[high]); 
    return (i + 1); 
} 

//Module Name: quickSort
//Author: https://www.geeksforgeeks.org/quick-sort/
//Date: March 12th, 2019
//Purpose: The main function used to quicksort our array of polygons
//Parameters: arr[]: the array to be sorted, low: the lowest point of the sort, high: the max size of the array
void quickSort(struct polygon arr[], int low, int high){ 
    if (low < high) 
    { 
        /* pi is partitioning index, arr[p] is now 
           at right place */
        int pi = partition(arr, low, high); 
  
        // Separately sort elements before 
        // partition and after partition 
        quickSort(arr, low, pi - 1); 
        quickSort(arr, pi + 1, high); 
    } 
} 

//Module Name: generateShapePolys
//Author: Zachary Kucera
//Date: March 12th, 2019
//Purpose: Uses the four points it recieves to create a polygon, and then calculates all the different light intensities, the normal, the distance from the camera, and sets the colour of the polygon.
//Parameters P0,P1,P2,P3: The four points that form the polygon, L: The light source matrix, E: The Camera position, C: the camera matrix, R,G,B: the red,green and blue components of the polygon's color
//Returns the fully constructed polygon
struct polygon generateShapePolys(dmatrix_t P0,dmatrix_t P1,dmatrix_t P2,dmatrix_t P3, dmatrix_t L, dmatrix_t E, dmatrix_t C, int R, int G, int B){
    dmatrix_t v1;
    dmatrix_t v2;
    
    v1 = *dmat_sub(from_homogeneous(&P1), from_homogeneous(&P0));//
    v2 = *dmat_sub(from_homogeneous(&P2), from_homogeneous(&P1));//Caclculate two vectors with the points

   
    p.world_points[0] = P0; //
    p.world_points[1] = P1; //  Create the polygon with the 4 world_points found above.
    p.world_points[2] = P2; //  
    p.world_points[3] = P3; //

    p.centroid = *dmat_scalar_mult(dmat_add(&P0, dmat_add(&P1, dmat_add(&P2, &P3))),0.25); //Find the centroid of the polygon

    p.normal = *dmat_normalize(dcross_product(&v1,&v2));//Find the normal of the poly and normalize it
    p.s = *dmat_normalize(dmat_sub(&L, &p.centroid));//Find the s vector and normalize it
    p.Id = (Ls*Pd)*max(0,((ddot_product(from_homogeneous(&p.s),&p.normal))/(dmat_norm(&p.s) * dmat_norm(&p.normal))));//Calculate the intensity of diffuse light for the polygon
    
    p.r = *dmat_normalize(dmat_add(dmat_scalar_mult(from_homogeneous(&p.s),-1), dmat_scalar_mult(&p.normal,(2*ddot_product(from_homogeneous(&p.s),&p.normal)/pow(dmat_norm(&p.normal),2)))));//Calculate the r vector and normalize it
    p.v = *dmat_normalize(dmat_sub(&E, &p.centroid));//Calculate the v vector and normailize it
    
    p.Is = (Ls * Ps) * max(0,(ddot_product(&p.r,from_homogeneous(&p.v)))/(dmat_norm(&p.r) * dmat_norm(&p.v)));

    p.distanceFromCamera = dmat_norm(dmat_sub(from_homogeneous(&E), from_homogeneous(&p.centroid)));//Calculate the distance from the camera
    p.camera_points[0] = *perspective_projection(dmat_mult(&C, &P0));//
    p.camera_points[1] = *perspective_projection(dmat_mult(&C, &P1));//Convert each point from world coordinates to 2D screen coordinates
    p.camera_points[2] = *perspective_projection(dmat_mult(&C, &P2));//
    p.camera_points[3] = *perspective_projection(dmat_mult(&C, &P3));//

    p.RED = R;
    p.GREEN = G;
    p.BLUE = B;
    
    return p;        
}

//Module Name: generateSpherePoints
//Author: Zachary Kucera
//Date: March 12th, 2019
//Purpose: Uses the parametric equation of a cone to construct all of the polygons needed to draw a Sphere.
//Parameters L: The light source matrix, E: The Camera position, C: the camera matrix, polygons: the array of polygons, count: the current number of polygons in polygons
//Returns count, so we can use the updated count in the next function
int generateSpherePoints(dmatrix_t L, dmatrix_t E, dmatrix_t C, struct polygon *polygons, int count){
    dmatrix_t P0;
    dmatrix_t P1;
    dmatrix_t P2;
    dmatrix_t P3;

    dmat_alloc(&P0,4,1) ;
    dmat_alloc(&P1,4,1) ;
    dmat_alloc(&P2,4,1) ;
    dmat_alloc(&P3,4,1) ;   
    
    static float dt = M_PI / 230;
    for (float u = 0.0; u <= M_PI; u += dt){ //Iterate u from 0 to PI
        for (float v = 0.0; v<= 2.0*M_PI + dt; v += dt){ //Iterate from v to 2PI
            P0.m[1][1] = (sin(u) * cos(v));  //X 3D coordinate of the sphere
            P0.m[2][1] = (sin(u) * sin(v));  //Y 3D coordinate of the sphere
            P0.m[3][1] = cos(u);             //Z 3D coordinate of the sphere
            P0.m[4][1] = 1.0;                //1 becuase parametric

            P1.m[1][1] = (sin(u + dt) * cos(v));  //X 3D coordinate of the sphere
            P1.m[2][1] = (sin(u + dt) * sin(v));  //Y 3D coordinate of the sphere
            P1.m[3][1] = cos(u + dt);             //Z 3D coordinate of the sphere
            P1.m[4][1] = 1.0;                //1 becuase parametric
            
            P2.m[1][1] = (sin(u + dt) * cos(v + dt));  //X 3D coordinate of the sphere
            P2.m[2][1] = (sin(u + dt) * sin(v+ dt));  //Y 3D coordinate of the sphere
            P2.m[3][1] = cos(u + dt);             //Z 3D coordinate of the sphere
            P2.m[4][1] = 1.0;                //1 becuase parametric

            P3.m[1][1] = (sin(u) * cos(v + dt));  //X 3D coordinate of the sphere
            P3.m[2][1] = (sin(u) * sin(v + dt));  //Y 3D coordinate of the sphere
            P3.m[3][1] = cos(u);             //Z 3D coordinate of the sphere
            P3.m[4][1] = 1.0;                //1 becuase parametric
        
            polygons[count] = generateShapePolys(P0,P1,P2,P3,L,E,C, 0,255,0);
            count ++;
        }
    }
    return count;
}

//Module Name: generateTorusPoints
//Author: Zachary Kucera
//Date: March 12th, 2019
//Purpose: Uses the parametric equation of a cone to construct all of the polygons needed to draw a torus.
//Parameters L: The light source matrix, E: The Camera position, C: the camera matrix, polygons: the array of polygons, count: the current number of polygons in polygons
//Returns count, so we can use the updated count in the next function
int generateTorusPoints(dmatrix_t L, dmatrix_t E, dmatrix_t C, struct polygon *polygons, int count){
    dmatrix_t P0;
    dmatrix_t P1;
    dmatrix_t P2;
    dmatrix_t P3;

    dmat_alloc(&P0,4,1) ;
    dmat_alloc(&P1,4,1) ;
    dmat_alloc(&P2,4,1) ;
    dmat_alloc(&P3,4,1) ;   

    float c = 3;    //How big the hole in the middle of the torus is
    float a = 0.7;  //Radius of the tube
    static float dt = M_PI / 195;

    for (float u = 0.0; u <= 2.0*M_PI; u += dt){            //Iterate both u and v to 2PI
        for (float v = 0.0; v<= 2.0*M_PI + dt; v += dt){    //
            P0.m[1][1] = ((c + (a * cos(v))) * cos(u));//
            P0.m[2][1] = ((c + (a * cos(v))) * sin(u));//Parametric Equation for a torus 
            P0.m[3][1] = (a * sin(v));       
            P0.m[4][1] = 1.0;

            P1.m[1][1] = ((c + (a * cos(v))) * cos(u + dt));//
            P1.m[2][1] = ((c + (a * cos(v))) * sin(u + dt));//Parametric Equation for a torus 
            P1.m[3][1] = (a * sin(v));       
            P1.m[4][1] = 1.0;  

            P2.m[1][1] = ((c + (a * cos(v + dt))) * cos(u + dt));//
            P2.m[2][1] = ((c + (a * cos(v + dt))) * sin(u + dt));//Parametric Equation for a torus 
            P2.m[3][1] = (a * sin(v + dt));       
            P2.m[4][1] = 1.0;  

            P3.m[1][1] = ((c + (a * cos(v + dt))) * cos(u));//
            P3.m[2][1] = ((c + (a * cos(v + dt))) * sin(u));//Parametric Equation for a torus 
            P3.m[3][1] = (a * sin(v + dt));       
            P3.m[4][1] = 1.0; 

            polygons[count] = generateShapePolys(P0,P1,P2,P3,L,E,C,255,0,0);
            count ++;
        }
    }
    return count;
}

//Module Name: generateConePoints
//Author: Zachary Kucera
//Date: March 12th, 2019
//Purpose: Uses the parametric equation of a cone to construct all of the polygons needed to draw a cone.
//Parameters L: The light source matrix, E: The Camera position, C: the camera matrix, polygons: the array of polygons, count: the current number of polygons in polygons
//Returns count, so we can use the updated count in the next function
int generateConePoints(dmatrix_t L, dmatrix_t E, dmatrix_t C, struct polygon *polygons, int count){
    dmatrix_t P0;
    dmatrix_t P1;
    dmatrix_t P2;
    dmatrix_t P3;

    dmat_alloc(&P0,4,1) ;
    dmat_alloc(&P1,4,1) ;
    dmat_alloc(&P2,4,1) ;
    dmat_alloc(&P3,4,1) ;   

    static float dt = M_PI / 100;
    static double dv = 0.002;
    for (float v = 0.0; v<= 1.0; v += dv){            //Iterate both u and v to 2PI
        for (float u = 0.0; u <= 2.0*M_PI; u += dt){    //
            P0.m[1][1] = (v-1.0) * cos(u) ;
            P0.m[2][1] = (v-1.0) * sin(u) ;//Parametric Equation for a cone 
            P0.m[3][1] = v + 1.1;       
            P0.m[4][1] = 1.0;

            P1.m[1][1] = (v-1.0) * cos(u + dt);
            P1.m[2][1] = (v-1.0) * sin(u + dt) ;//Parametric Equation for a cone 
            P1.m[3][1] = v + 1.1;       
            P1.m[4][1] = 1.0;

            P2.m[1][1] = (v -1.0 + dv) * cos(u + dt);
            P2.m[2][1] = (v -1.0 + dv) * sin(u + dt);//Parametric Equation for a cone 
            P2.m[3][1] = v + dv + 1.1;       
            P2.m[4][1] = 1.0;

            P3.m[1][1] = (v -1.0 + dv) * cos(u);
            P3.m[2][1] = (v -1.0 + dv) * sin(u);//Parametric Equation for a cone 
            P3.m[3][1] = v + dv + 1.1;       
            P3.m[4][1] = 1.0;

            polygons[count] = generateShapePolys(P0,P1,P2,P3,L,E,C,0,255,255);
            count ++;
        }
    }
    return count;
}



//Module Name: Draw
//Author: Zachary Kucera
//Date: March 12th, 2019
//Purpose: The primary function called to draw the shapes. Defines our camera, sorts polygons, and calls the XFil function to fill the polygons
void draw() {
    hdc = GetDC(hwnd);

    dmatrix_t E ; /* The centre of projection for the camera */
    
    dmat_alloc(&E,4,1) ;
    
    E.m[1][1] = Ex ;
    E.m[2][1] = Ey ;
    E.m[3][1] = Ez ;
    E.m[4][1] = 1.0 ;
    
    dmatrix_t G ; /* Point gazed at by camera */
    
    dmat_alloc(&G,4,1) ;
    
    G.m[1][1] = Gx ;
    G.m[2][1] = Gy ;
    G.m[3][1] = Gz ;
    G.m[4][1] = 1.0 ;

    dmatrix_t C ; /* The camera matrix */

    dmat_alloc(&C,4,4) ;
    C = *build_camera_matrix(&E,&G) ;

    /* The light source matrix*/
    dmatrix_t L;
    dmat_alloc(&L,4,1);
    L.m[1][1] = Lx ;
    L.m[2][1] = Ly ;
    L.m[3][1] = Lz ;
    L.m[4][1] = 1.0;  

    static struct polygon polygons[359412]; //The array that contains all of the polygons for the various shapes
    int count = 0;//Keeps track of how many polygons we have

    count = generateSpherePoints(L,E,C, polygons, count);// This adds the sphere polys to the array, returns count so we know how many polys we have
    printf("\nPOST SPHERE: %d", count);

    count = generateTorusPoints(L,E,C, polygons, count);// This adds the torus polys to the array, returns count so we know how many polys we have
    printf("\nPOST TORUS: %d", count);

    count = generateConePoints(L,E,C, polygons, count);// This adds the sphere cone to the array, returns count so we know how many polys we have
    printf("\nPOST CONE: %d", count);
    
    float I;//The total light intensity for any given polygon
    quickSort(polygons,0,count);//Sort the polygons by distance from the camera

    for(int i = 0; i < count; i++){
        I = polygons[i].Id + Ia * Pa + polygons[i].Is;//Add up the three different types of light to get the total light intensity.
        XFillConvexPolygon(hdc, RGB((int)polygons[i].RED* I,(int)polygons[i].GREEN*I,(int)polygons[i].BLUE* I), polygons[i].camera_points, 4); //fill the polys, using thier I value to determine the intensity of the colour
    }
    
}
//Module Name: WndProc
//Author: http://www.winprog.org/tutorial/simple_window.html added upon by Zachary Kucera
//Date: Jan 15th, 2019
//Purpose: Handels the different messages sent by the window. This is what begins the painting process
LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
    hdc = GetDC(hwnd);

    switch(msg)
    {
		case WM_PAINT: //When a WM_PAINT message is recieved, begin paint, draw the shape, and end paint.
            hdc = BeginPaint(hwnd, &ps);
            draw();
            EndPaint(hwnd, &ps);
        break;

        case WM_CHAR://When a 'q' is pressed, close the window.
            if (wParam == 113) DestroyWindow(hwnd);
        break;

        case WM_CLOSE:// When a close message is recieved, close the window
            DestroyWindow(hwnd);
        break;
        case WM_DESTROY:
            PostQuitMessage(0);
        break;
        default:
            return DefWindowProc(hwnd, msg, wParam, lParam);
    }
    return 0;
}


//Module Name: WinMain
//Author: http://www.winprog.org/tutorial/simple_window.html 
//Date: Jan 15th, 2019
//Purpose: Main window function, opens the window and sends messages.
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
    LPSTR lpCmdLine, int nCmdShow)
{
    //Step 1: Registering the Window Class
    wc.cbSize        = sizeof(WNDCLASSEX);
    wc.style         = 0;
    wc.lpfnWndProc   = WndProc;
    wc.cbClsExtra    = 0;
    wc.cbWndExtra    = 0;
    wc.hInstance     = hInstance;
    wc.hIcon         = LoadIcon(NULL, IDI_APPLICATION);
    wc.hCursor       = LoadCursor(NULL, IDC_ARROW);
    wc.hbrBackground = (HBRUSH)(COLOR_WINDOW+1);
    wc.lpszMenuName  = NULL;
    wc.lpszClassName = g_szClassName;
    wc.hIconSm       = LoadIcon(NULL, IDI_APPLICATION);

    if(!RegisterClassEx(&wc)) //If the window is not properly created
    {
        MessageBox(NULL, "Window Registration Failed!", "Error!",
            MB_ICONEXCLAMATION | MB_OK);
        return 0;
    }

    // Create the window
    hwnd = CreateWindowEx(
        WS_EX_CLIENTEDGE,
        g_szClassName,
        "Computer Graphics Assignment 3",
        WS_OVERLAPPEDWINDOW,
        CW_USEDEFAULT, CW_USEDEFAULT, 512, 512, // Set the size of the window to 512 x 512
        NULL, NULL, hInstance, NULL);

    if(hwnd == NULL)
    {
        MessageBox(NULL, "Window Creation Failed!", "Error!",
            MB_ICONEXCLAMATION | MB_OK);
        return 0;
    }

    ShowWindow(hwnd, nCmdShow);//Display the window
    UpdateWindow(hwnd);

    // This is the message loop. This dispatches the messages to WndProc, where the messages are handled.
    while(GetMessage(&Msg, NULL, 0, 0) > 0)
    {
        TranslateMessage(&Msg);
        DispatchMessage(&Msg);
    }
    return Msg.wParam;
}

