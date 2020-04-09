// lab2Parallels.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <omp.h>
#include <Windows.h>
using namespace std;

 
void task1() {
#pragma omp parallel num_threads(4)
    {
        printf("Hello World!\n");
    }
}
void task2() {
#pragma omp parallel num_threads(3)
    {
        int k = omp_get_thread_num();
        if (k % 2 == 0) {
            printf("I am %d thread from %d threads!\n",
                k,
                omp_get_num_threads()
            );
        }
    }
}
void task3() {
    int rank;
    omp_set_num_threads(3);
#pragma omp parallel private(rank)
    {
        rank = omp_get_thread_num();
        Sleep(100);
        printf("I am %d thread.\n", rank);
    }
}
void task4() {
    int sum = 0;
    int n = 4;
#pragma omp parallel num_threads(2) reduction(+:sum)
    {
        if (omp_get_thread_num() == 0) {
            for (int i = 1; i <= n / 2; i++)
                sum += i;
        }
        else {
            for (int i = n / 2; i <= n; i++)
                sum += i;
        }
    }
    printf("sum = %d\n", sum);
}
void task5() {
    omp_set_num_threads(2);
    int N = 4;
    int sum = 0;
#pragma omp parallel
    {
#pragma omp for reduction(+:sum)
        for (int i = 1; i <= N; i++) {
            sum += i;
        }
    }
    printf("sum = %d\n", sum);
}
void task6() {
    omp_set_num_threads(4);
    int N = 10;
    int sum = 0;
#pragma omp parallel
    {
#pragma omp for reduction(+:sum) schedule(guided, 2)
        for (int i = 1; i <= N; i++) {
            sum += i;
            printf("%d:calculation of the iteration number %d.\n", omp_get_thread_num(), i);
        }
    }
    printf("sum = %d\n", sum);
}
void task7() {
    int n = 1000000000;
    double pi = 0, xi;
#pragma omp parallel
    {
#pragma omp for reduction(+:pi)
        for (int i = 0; i < n; i++) {
            xi = (i + 0.5)/ n;
            pi += 4 /( 1 + xi * xi);
        }
    }
    pi /= n;
    std::cout << pi;
}
void task8() {
    int N;
    cout << "\nN:\n";
    cin >> N;
    double **A = new double*[N];
    double **B = new double* [N];
    double **C = new double* [N];
    for (int i = 0; i < N; i++) {
        A[i] = new double[N];
        B[i] = new double[N];
        C[i] = new double[N];
    }
    cout << "Matrix A:\n";
    for (int i = 0; i < N; i++) 
        for (int j = 0; j < N; j++)
            cin >> A[i][j];
    cout << "Matrix B:\n";
    for (int i = 0; i < N; i++) 
        for (int j = 0; j < N; j++)
            cin >> B[i][j];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C[i][j] = 0;
#pragma omp parallel
            {
#pragma omp for
                for (int k = 0; k < N; k++)
                    C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    cout << "Matrix C:\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++)
            cout << C[i][j] << " ";
        cout << endl;
    }
}
void task9() {
    int k, ki;
    cout << "k:\n";
    cin >> k;
    if (k > 4)
        k = 4;

#pragma omp parallel
{   
    for (int ki = 1; ki <= k; ki++) {
#pragma omp sections
        {
            printf("%d: parallel region\n", ki);
#pragma omp section
            {
                printf("%d : came in section 1\n", ki);
            }
#pragma omp section
            {
                printf("%d : came in section 2\n", ki);
            }
#pragma omp section
            {
                printf("%d : came in section 3\n", ki);
            }
        }
    }
}
}
void task10() {
    omp_set_num_threads(2);
    int N = 4;
    int sum = 0;
#pragma omp parallel
    {
#pragma omp for
        for (int i = 1; i <= N; i++) {
#pragma omp atomic
            sum += i;
        }
    }
    printf("sum = %d\n", sum);
}
void task11() {
    int n = 10000000;
    double pi = 0;
    omp_set_num_threads(4);
#pragma omp parallel
    {
#pragma omp for schedule(static)
        for (int i = 0; i < n; i++) {
#pragma omp critical(omp)
            {
                pi += 4 / (1 + ((i + 0.5) / n) *((i + 0.5) / n));
            }
        }
    }
    pi /= n;
    std::cout << pi;
}

int main()
{
   /* task1();
    task2();
    task3();
    task4();
    task5();*/
    //task6();
    //double start = omp_get_wtime();
    task7();
    //double end = omp_get_wtime();
    //double wtick = omp_get_wtick();
    /*printf("Time execution = %.16g\n", end - start);
    printf("wtick = %.16g\n", wtick);
    */task8();
    task9();
    task10();
   // start = omp_get_wtime();
    task11();
   // end = omp_get_wtime();
    //wtick = omp_get_wtick();
   // printf("Time execution = %.16g\n", end - start);
    //printf("wtick = %.16g\n", wtick);*/
    return 0;
}
   
