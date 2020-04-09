// Lab3Parallels.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "mpi.h"
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <iostream>
using namespace std;

void task1() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    printf("I am %d process from %d processes!\n", rank, size);
}

void task2() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        int size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        printf("%d processes!\n", size);
    }
    else if(rank % 2 == 0) {
        printf("I am %d: SECOND!", rank);
    }
    else{
        printf("I am %d: FIRST!", rank);
    }
}

void task3() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0){
        int buf = 10;
        MPI_Send(&buf, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    }
    else {
        int buf;
        MPI_Status status;
        MPI_Recv(&buf, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        printf("this process %d, i get message - %d", rank, buf);
    }
}

void task4() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int buf;
    MPI_Status status;
    if (rank == 0) {
        buf = 0;
        MPI_Send(&buf, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        MPI_Recv(&buf, 1, MPI_INT, size - 1, 0, MPI_COMM_WORLD, &status);
        printf("this process %d, i get message - %d", rank, buf);
    }
    else {
        MPI_Recv(&buf, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &status);
        printf("this process %d, i get message - %d", rank, buf);
        buf++;
        if (rank == size - 1)
            MPI_Send(&buf, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        else
            MPI_Send(&buf, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    }
}

void task5() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int buf;
    MPI_Status status;
    if (rank == 0) {
        for (int i = 1; i <= size - 1; i++) {
            MPI_Recv(&buf, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
            printf("this process %d, i get message - %d\n", rank, buf);
       }
    }
    else {
        buf = rank;
        MPI_Send(&buf, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
}

void task6() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int buf;
    MPI_Request reqSend, reqRecv;
    MPI_Status statSend, statRecv;
    if (rank == 0) {
        buf = 200;
        MPI_Isend(&buf, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &reqSend);
        MPI_Wait(&reqSend, &statSend);
    }
    else {
        MPI_Irecv(&buf, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &reqRecv);
        MPI_Wait(&reqRecv, &statRecv);
        printf("receive message %d\n", buf);
    }
}

void task7()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int buf;
    MPI_Request reqSend, reqRecv;
    MPI_Status statSend, statRecv;
    if (rank == size - 1) {
        buf = rank;
        MPI_Isend(&buf, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &reqSend);
        MPI_Irecv(&buf, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &reqRecv);
        MPI_Wait(&reqSend, &statSend);
        MPI_Wait(&reqRecv, &statRecv);
        printf("[%i] receive message %d\n", rank, buf);
    }
    else if (rank == 0)
    {
        buf = rank;
        MPI_Isend(&buf, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD, &reqSend);
        MPI_Irecv(&buf, 1, MPI_INT, size - 1, 0, MPI_COMM_WORLD, &reqRecv);
        MPI_Wait(&reqSend, &statSend);
        MPI_Wait(&reqRecv, &statRecv);
        printf("[%i] receive message %d\n", rank, buf);
    }
    else {
        buf = rank;
        MPI_Isend(&buf, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD, &reqSend);
        MPI_Irecv(&buf, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &reqRecv);
        MPI_Wait(&reqSend, &statSend);
        MPI_Wait(&reqRecv, &statRecv);
        printf("[%i] receive message %d\n", rank, buf);
    }
}

void task8() {
    int rank;
    MPI_Request req[4];
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int buf[3];
    for (int i = 0; i <= size - 1; i++) {
        if (i != rank) {
            MPI_Isend(&rank, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &req[i]);
        }
        else {
            req[i] = MPI_REQUEST_NULL;
        }
    }
    MPI_Waitall(4, req, MPI_STATUS_IGNORE);
    for (int i = 0, j = 0; i <= size - 1; i++) {
        if (i != rank) {
            MPI_Irecv(&buf[j++], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &req[i]);
        }
        else {
            req[i] = MPI_REQUEST_NULL;
        }
    }
    MPI_Waitall(4, req, MPI_STATUS_IGNORE);
    for(int i = 0; i<3; i++)
        printf("this process %d, i get message - %d\n", rank, buf[i]);
}

void task9() {
    int rank;
    MPI_Request req[4];
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    char buf[100];
    int n;
    if (rank == 0) {
        scanf_s("%i", &n);
        strcpy_s(buf, "rtvabc");

    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(buf, n, MPI_CHAR, 0, MPI_COMM_WORLD);
    char a = 'a';
   
    printf("rtvabc\n");

    for (int i = rank; i < 26; i += size) {
        int count = 0;
        int num = (int)a + i;
        for (int j = 0; j < n; j++) {
            if ((char)num == buf[j])
                count++;
        }
        if(count>0)
            printf("%c, count = %d\n", char(num), count);
    }


}
void task10()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    long double drob, drobSum = 0, Result = 0;

    int n = 1000000;

    for (int i = rank; i <= n; i += size)
    {
        drob = 4 / (1.0 + pow((i + 0.5) * (1.0 / n), 2.0));
        drobSum += drob;
    }

    MPI_Reduce(&drobSum, &Result, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        long double x = Result * (1.0 / n);
        printf("Pi is %lf", x);
    }

}

/*void task11()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int n;
    int elem_per_proc;
    int** a, ** b, ** c, ** sendBuf;

    a = new int* [1];
    b = new int* [1];
    c = new int* [1];
    sendBuf = new int* [1];


    if (rank == 0)
    {
        scanf_s("%d", &n);
        int buf;
        flushall();

        a = new int* [n];
        for (int i = 0; i < n; i++) {
            a[i] = new int[n];
        }

        b = new int* [n];
        for (int i = 0; i < n; i++) {
            b[i] = new int[n];
        }

        c = new int* [n];
        for (int i = 0; i < n; i++) {
            c[i] = new int[n];
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                printf("Input value for %d-%d elem in a", i, j);
                scanf_s("%d", &a[i][j]);
                flushall();
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                printf("Input value for %d-%d elem in b", i, j);
                scanf_s("%d", &b[i][j]);
                flushall();
            }

            elem_per_proc = n / size;
            if (elem_per_proc == 0)elem_per_proc++;
            else if (n % size >= 1)
            {
                elem_per_proc++;
            }
        }
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    sendBuf = new int* [n];
    for (int i = 0; i < n; i++) {
        sendBuf[i] = new int[n];
    }


    MPI_Bcast(&elem_per_proc, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(b, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Scatter(a, elem_per_proc, MPI_INT, sendBuf, elem_per_proc, MPI_INT, 0, MPI_COMM_WORLD);

    int** buf_mas;
    buf_mas = new int* [elem_per_proc];
    for (int i = 0; i < elem_per_proc; i++) {
        c[i] = new int[n];
    }

    for (int i = 1; i <= elem_per_proc; i++)
    {
        for (int j = 0; j < n; j++)
            buf_mas[i][j] = sendBuf[i][j] * b[j][rank * i];
    }

    MPI_Gather(buf_mas, 1, MPI_INT, c, 1, MPI_INT, 0,
        MPI_COMM_WORLD);

    if (rank == 0)
    {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                printf("%d-%d value of C is %d", i, j, c[i][j]);

            }
        }
    }
}

void task12() {
    MPI_Group group, oldgroup;
    MPI_Comm newcomm;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int ranks[2] = { 0, 2 };
    MPI_Comm_group(MPI_COMM_WORLD, &oldgroup);
    MPI_Group_incl(oldgroup, 2, ranks, &group);
    MPI_Comm_create(MPI_COMM_WORLD, group, &newcomm);
    char message = '0';
    int newrank = 0;
    int newsize = 0;
    if (newcomm != MPI_COMM_NULL) {
        MPI_Comm_rank(newcomm, &newrank);
        MPI_Comm_size(newcomm, &newsize);
        if (newrank == 0) {
            message = 'A';
        }
        MPI_Bcast(&message, 1, MPI_CHAR, 0, newcomm);
    }
    printf("MPI_COMM_WORLD: %d from %d. New comm : %d from %d.Message = %c.", rank, size, newrank, newsize, message);
}*/


void KeyboardDataInit(double* m, double* v, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            cout << "a[" << i + 1 << "]" << '[' << j + 1 << "]=";
            cin >> m[i * size + j];
        }
        cout << "b[" << i + 1 << "]= ";
        cin >> v[i];
    }
}
void RandomDataInit(double* m, double* v, int size) {
    srand(time(0));
    for (int i = 0; i < size; i++) {
        double sum = 0;
        for (int j = 0; j < size; j++) {
            if (i != j) {
                m[i * size + j] = rand() % 18 - 9;
                sum += abs(m[i * size + j]);
            }
        }
        int sign;
        if (rand() > 16000)
            sign = 1;
        else
            sign = -1;
        m[i * size + i] = sign * sum * (rand() / 32767.0 + 1);
        v[i] = rand() % 18 - 9;
    }
}

int main(int argc, char* argv[])
{
    double* pA;
    double* pB;
    double* pX;
    double* pProcRowsA;
    double* pProcRowsB;
    double* pProcTempX;
    int Size, RowNum, ProcNum, ProcRank;
    double Eps, Norm, MaxNorm;
    double Start, Finish, Dt;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    if (ProcRank == 0) {
        cout << "Input the dimension  ";
        cin >> Size;
    }
    MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    RowNum = Size / ProcNum;
    Eps = 0.001;
    pX = new double[Size];
    pA = 0;
    pB = 0;
    pProcRowsA = new double[RowNum * Size];
    pProcRowsB = new double[RowNum];
    pProcTempX = new double[RowNum];
    if (ProcRank == 0) {
        pA = new double[Size * Size];
        pB = new double[Size];
        RandomDataInit(pA, pB, Size);
        //KeyboardDataInit(pA, pB, Size);
        for (int i = 0; i < Size; i++) {
            pX[i] = 0;
        }
    }
    MPI_Bcast(pX, Size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(pA, RowNum * Size, MPI_DOUBLE, pProcRowsA, RowNum * Size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(pB, RowNum, MPI_DOUBLE, pProcRowsB, RowNum, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    Start = MPI_Wtime();
    do {
        for (int i = 0; i < RowNum; i++) {
            pProcTempX[i] = -pProcRowsB[i];
            for (int j = 0; j < Size; j++) {
                if (ProcRank * RowNum + i != j)
                    pProcTempX[i] += pProcRowsA[i * Size + j] * pX[j];
            }
            pProcTempX[i] /= -pProcRowsA[i * Size + i];
        }
        Norm = fabs(pX[ProcRank * RowNum] - pProcTempX[0]);
        for (int i = 0; i < RowNum; i++) {
            if (fabs(pX[ProcRank * RowNum + i] - pProcTempX[i]) > Norm)
                Norm = fabs(pX[ProcRank * RowNum + i] - pProcTempX[i]);
        }
        MPI_Reduce(&Norm, &MaxNorm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Bcast(&MaxNorm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Allgather(pProcTempX, RowNum, MPI_DOUBLE, pX, Size, MPI_DOUBLE, MPI_COMM_WORLD);
    } while (MaxNorm > Eps);
    Finish = MPI_Wtime();
    if (ProcRank == 0) {
        cout << endl << "time solution of system=" << (Finish - Start) << " sec." << endl;
        for( int i=0;i<Size;i++){
            cout.width(5);
            cout.precision(4);
            cout<<"x["<<i+1<<"]= "<<pX[i]<<endl;
        }
    }
    delete[]pA;
    delete[]pB;
    delete[]pX;
    delete[]pProcRowsA;
    delete[]pProcRowsB;
    delete[]pProcTempX;
    MPI_Finalize();
    /*MPI_Init(&argc, &argv);
    printf("2\n1 3\n4 8\n5 4\n3 0\n");
    printf("result:\n14 4\n44 16");


    MPI_Finalize();*/
    return 0;
}


