using System;
using System.Linq;
using System.Text;

class MatrixOperations
{
    static MatrixOperations()
    {
        // Встановлення UTF-8 кодування для підтримки українських символів та символу "λ"
        Console.OutputEncoding = Encoding.UTF8;
    }

    // Перевірка умов матриці: симетричність і визначники мінорів
    public static bool CheckMatrixConditions(double[,] A)
    {
        int n = A.GetLength(0);
        bool isSymmetric = IsSymmetric(A);
        bool minorsPositive = true;

        for (int i = 1; i < n; i++)
        {
            double[,] minor = GetMinor(A, i + 1);
            double detMinor = Determinant(minor);
            detMinor = Math.Round(detMinor, 15);
            string sign = detMinor > 0 ? ">" : "<=";
            Console.WriteLine($"Визначник мінору M{i + 1}x{i + 1} = {detMinor} ({sign} 0)");

            if (detMinor <= 0)
            {
                minorsPositive = false;
            }
        }

        return isSymmetric && minorsPositive;
    }

    // Метод степеневої ітерації для знаходження максимального власного значення
    public static double PowerIteration(double[,] A, double eps = 1e-4, int maxIterations = 100)
    {
        int n = A.GetLength(0);
        /* double[] x = Enumerable.Repeat(1.0, n).ToArray(); // початковий вектор */
        double[] x = new double[] { -10.0, 1.0, 10.0, 1.0 }; 
        /* double[] x = Enumerable.Repeat(-10.0, n).ToArray(); // початковий вектор */
        double lambda1 = 0;

        for (int i = 0; i < maxIterations; i++)
        {
            double[] xNext = MultiplyMatrixVector(A, x);

            if (i == 0)
            {
                lambda1 = xNext[0];
                Console.WriteLine($"Ітерація {i + 1} - λ1 = {lambda1}");
            }
            else
            {
                double normXNext = Math.Sqrt(xNext.Select(xi => xi * xi).Sum()); // евклідова норма вектора
                double[] xNorm = xNext.Select(xi => xi / normXNext).ToArray();

                double[] aXNorm = MultiplyMatrixVector(A, xNorm);
                double lambda1Next = aXNorm[0] / xNorm[0];
                Console.WriteLine($"Ітерація {i + 1} : λ1 = {lambda1Next:F6}".PadRight(45) + $"|λ_k+1 - λ_k| = {Math.Abs(lambda1Next - lambda1):F6}");

                if (Math.Abs(lambda1Next - lambda1) <= eps)
                {
                    Console.WriteLine($"Ітераційний процес зупинено. λ1 = {lambda1Next:F6} (ε = {eps})");
                    break;
                }

                lambda1 = lambda1Next;
                x = xNorm;
            }
        }

        return lambda1;
    }

    public static void Main(string[] args)
    {
        double eps = 1e-4;

        double[,] A = {
             { 10.0, 2.0, 1.0, 2.0 },
             { 2.0, 10.0, 2.0, 1.0 },
             { 1.0, 2.0, 10.0, 1.0 },
             { 2.0, 1.0, 1.0, 10.0 }
         };

        Console.WriteLine("Вхідні дані:");
        PrintMatrix(A);
        Console.WriteLine($"Точність ε = {eps}\n");

        bool matrixConditions = CheckMatrixConditions(A);

        if (matrixConditions)
        {
            Console.WriteLine("Матриця A є симетричною і всі визначники мінорів додатні.");
            Console.WriteLine("Отже, можна знайти мінімальне власне значення як λmin(A) = ||A||∞ − λmax(B). \n");
        }
        else
        {
            Console.WriteLine("Не виконуються умови для знаходження мінімального власного значення.");
            return;
        }

        double normAInf = InfinityNorm(A);
        Console.WriteLine($"Норма матриці за ||A||inf нормою: {normAInf}");

        Console.WriteLine("\n--Розрахунок λmax(A)--");
        double maxEigenvalueA = PowerIteration(A, eps);
        Console.WriteLine($"\nМаксимальне власне значення λmax(A) = {maxEigenvalueA}");

        if (matrixConditions)
        {
            Console.WriteLine("\n--Розрахунок λmax(B)--");
            double[,] B = SubtractMatrixFromIdentity(A, maxEigenvalueA);
            double maxEigenvalueB = PowerIteration(B, eps);
            Console.WriteLine($"\nМаксимальне власне значення λmax(B) = {maxEigenvalueB}");

            double lambdaMinA = normAInf - maxEigenvalueB;
            Console.WriteLine($"Мінімальне власне значення λmin(A) = ||A||∞ − λmax(B) = {lambdaMinA}");
        }
    }

    // Допоміжні функції для обчислення мінору, перевірки симетричності, визначення та норм
    public static bool IsSymmetric(double[,] A)
    {
        int n = A.GetLength(0);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (A[i, j] != A[j, i]) return false;
        return true;
    }

    public static double Determinant(double[,] matrix)
    {
        if (matrix.GetLength(0) == 2)
            return matrix[0, 0] * matrix[1, 1] - matrix[0, 1] * matrix[1, 0];

        double det = 0;
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            double[,] subMatrix = GetSubMatrix(matrix, i);
            det += (i % 2 == 0 ? 1 : -1) * matrix[0, i] * Determinant(subMatrix);
        }
        return det;
    }

    public static double[,] GetMinor(double[,] A, int size)
    {
        double[,] minor = new double[size, size];
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                minor[i, j] = A[i, j];
        return minor;
    }

    public static double[,] GetSubMatrix(double[,] matrix, int excludeCol) 
    {
        int size = matrix.GetLength(0) - 1;
        double[,] result = new double[size, size];
        for (int i = 1; i < matrix.GetLength(0); i++)
            for (int j = 0, col = 0; j < matrix.GetLength(1); j++)
                if (j != excludeCol) result[i - 1, col++] = matrix[i, j];
        return result;
    }

    public static double InfinityNorm(double[,] A)
    {
        double maxSum = 0;
        for (int i = 0; i < A.GetLength(0); i++)
        {
            double rowSum = 0;
            for (int j = 0; j < A.GetLength(1); j++)
                rowSum += Math.Abs(A[i, j]);
            maxSum = Math.Max(maxSum, rowSum);
        }
        return maxSum;
    }

    public static double[] MultiplyMatrixVector(double[,] A, double[] x)
    {
        int n = A.GetLength(0);
        double[] result = new double[n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                result[i] += A[i, j] * x[j];
        return result;
    }

    public static double[,] SubtractMatrixFromIdentity(double[,] A, double scalar)
    {
        int n = A.GetLength(0);
        double[,] result = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                result[i, j] = (i == j ? scalar : 0) - A[i, j];
        return result;
    }

    public static void PrintMatrix(double[,] A)
    {
        int n = A.GetLength(0);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                Console.Write($"{A[i, j],7}");
            Console.WriteLine();
        }
    }
}
