//
// Created by lzz20 on 2022/9/26.
//
#include <iostream>
#include <omp.h>
#include <chrono>
#include <vector>

//ex01 a



long Fibonacci_sequential_recursive(int N)
{
    long a, b;

    if (N < 2)
    {
        return N;
    }

        a = Fibonacci_sequential_recursive(N - 1);

        b = Fibonacci_sequential_recursive(N - 2);

    return a + b;

}

long Fibonacci_parallel_recursive(int N)
{
    long a, b;

    if (N < 2)
        {
        return N;
        }

#pragma omp task shared(a) firstprivate(N)
    {
        a = Fibonacci_parallel_recursive(N - 1);
    }
#pragma omp task shared(b) firstprivate(N)
    {
        b = Fibonacci_parallel_recursive(N - 2);
    }
#pragma omp taskwait
    return a + b;
}

long Fibonacci_sequential_iterative(int N)
{
    std::vector<int>res={0,1};
    if (N<2)
    {
        return res[N];
    }

    int i=1;

    long fib1=0;
    long fib2=1;
    long fib=0;
    while (i<N) {

        fib=fib1+fib2;
        fib1=fib2;
        fib2=fib;
        i++;
    }

    return fib;
}

 long Fibonacci_parallel_iterative(int N)
 {
     std::vector<int>res={0,1};
     if (N<2)
     {
         return res[N];
     }

     int i=1;
     long fib1=0;
     long fib2=1;
     long fib=0;
#pragma omp parallel
     {
#pragma omp single nowait
         {
             while (i < N) {
/*#pragma omp task firstprivate(i) shared(fib,fib1,fib2) depend(in:fib1) depend(in:fib2) depend(out:fib)
                 {
                     fib = fib1 + fib2;

                 }
#pragma omp task firstprivate(i) shared(fib1,fib2) depend(in:fib2) depend(out:fib1)
                 {
                     fib1 = fib2;
                 }
#pragma omp task firstprivate(i) shared(fib,fib2) depend(in:fib) depend(out:fib2)
                 {
                     fib2 = fib;
                 }*/
#pragma omp task /*firstprivate(i)*/ shared(fib,fib1,fib2) depend(in:fib1) depend(in:fib2) depend(out:fib)
                 {
                     fib = fib1 + fib2;
                     fib1 = fib2;
                     fib2 = fib;
                 }
#pragma omp atomic
                 i++;
             }
         }
     }

     return fib;
 }


int main() {
    int thread_number = 8;

    omp_set_num_threads(thread_number);

    int n = 10;

    auto start = std::chrono::high_resolution_clock::now();
//sequential recursive execute
    std::cout << "fib " << n << " value = " << Fibonacci_sequential_recursive(n) << std::endl;
    std::chrono::duration<double> temps_seq_rec = std::chrono::high_resolution_clock::now() - start;
    std::cout << "Recursive Sequential Time: " << temps_seq_rec.count() << "s\n";

    start = std::chrono::high_resolution_clock::now();
    //parallel recursive execute
#pragma omp parallel shared(n)
    {
#pragma omp single nowait
        std::cout << "fib " << n << " value = " << Fibonacci_parallel_recursive(n) << std::endl;
    }
    std::chrono::duration<double> temps_para_rec = std::chrono::high_resolution_clock::now() - start;
    std::cout << "Recursive Parallel Time: " << temps_para_rec.count() << "s\n";

    std::cout << "Speed up recursive = " << temps_seq_rec.count()/temps_para_rec.count() << "\n";

    start = std::chrono::high_resolution_clock::now();
//sequential iterative execute
    std::cout << "fib " << n << " value = " << Fibonacci_sequential_iterative(n) << std::endl;
    std::chrono::duration<double> temps_seq_ite = std::chrono::high_resolution_clock::now() - start;
    std::cout << "Iterative sequential Time: " << temps_seq_ite.count() << "s\n";

    start = std::chrono::high_resolution_clock::now();
//parallel iterative execute
    std::cout << "fib " << n << " value = " << Fibonacci_parallel_iterative(n) << std::endl;
    std::chrono::duration<double> temps_para_ite = std::chrono::high_resolution_clock::now() - start;
    std::cout << "Iterative parallel Time: " << temps_para_ite.count() << "s\n";

    std::cout << "Speed up iterative = " << temps_seq_ite.count()/temps_para_ite.count() << "\n";
    std::cout << "Speed up compare with sequential recursive = " << temps_seq_rec.count()/temps_para_ite.count() << "\n";
























    return 0;
}