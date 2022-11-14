#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<cmath>
#define random(a,b) (rand()%(b-a+1)+a)
using namespace std;
class Matrix {
	public:
		long long int** m;//矩阵 
		long long int size;//阶数 
		Matrix() {
			m = NULL;
			size = 0;
		}
		Matrix(long long int size) {
			this->size = size;
			m = new long long int*[size];
			for (long long int i = 0; i < size; i++) {
				m[i] = new long long int[size];
				for (long long int j = 0; j < size; j++) {
					m[i][j] = 0;
				}
			}
		}
		//矩阵打印 
		void Matrix_Print(long long int** MatrixA, long long int size);
		//矩阵切分 
		Matrix Matrix_Split(long long int** MatrixA, long long int x, long long int y, long long int size);
		//判断矩阵阶数是否2的幂 
		bool IsPower(long long int size);
		//矩阵加法 
		void Matrix_Sum(long long int**  MatrixA, long long int** MatrixB, long long int** MatrixResult, long long int size);
		//矩阵减法 
		void Matrix_Sub(long long int**  MatrixA, long long int** MatrixB, long long int** MatrixResult, long long int size);
		//矩阵传统乘法 
		void Matrix_Mul(long long int**  MatrixA, long long int** MatrixB, long long int** MatrixResult, long long int size);
		//矩阵Strassen乘法 
		void StrassenMul(long long int**  MatrixA, long long int** MatrixB, long long int** MatrixResult, long long int size);
		//偶数非2的幂次阶矩阵乘法 
		void EvenNumberMul(long long int** matrixA,long long int** matrixB,long long int** MatrixResult,long long int size);
		//奇数阶矩阵乘法 
		void OddNumberMul(long long int** matrixA,long long int** matrixB,long long int** MatrixResult,long long int size);
		//任意阶矩阵乘法 
		void Matrix_Mul_AllType(long long int** matrixA,long long int** matrixB,long long int** matrixC,long long int size);
		//矩阵赋值 
		void MakeMatrix(long long int**  MatrixA, long long int size);
		//将矩阵各个元素求和，如果两种算法得出的矩阵的和相等则认为算法正确 
		long long int GetMatrixSum(long long int** Matrix, long long int size);
		
};

//矩阵打印 
void Matrix::Matrix_Print(long long int** MatrixA, long long int size) {
	for (long long int i = 0; i < size; i++) {
		for (long long int j = 0; j < size; j++)
			cout<<MatrixA[i][j]<<" ";
		cout<<endl;
	}
	cout<<endl;
}


//判断阶数是否为2的幂次
bool Matrix::IsPower(long long int size) {
	if(size<1) return false;
	long long int m = size&(size-1);
	return m == 0;
}

void Matrix::Matrix_Sum(long long int**  MatrixA, long long int** MatrixB, long long int** MatrixResult, long long int size) {
	for(long long int i = 0; i < size; i++) {
		for(long long int j = 0; j < size; j++) {
			MatrixResult[i][j] = MatrixA[i][j] + MatrixB[i][j];
		}
	}
}


void Matrix::Matrix_Sub(long long int**  MatrixA, long long int** MatrixB, long long int** MatrixResult, long long int size) {
	for(long long int i = 0; i < size; i++) {
		for(long long int j = 0; j < size; j++) {
			MatrixResult[i][j] = MatrixA[i][j] - MatrixB[i][j];
		}
	}
}

//传统矩阵乘法 
void Matrix::Matrix_Mul(long long int**  MatrixA, long long int** MatrixB, long long int** MatrixResult, long long int size) {
	for(long long int i = 0; i < size; i++) {
		for(long long int j = 0; j < size; j++) {
			MatrixResult[i][j] = 0;
			for(long long int k = 0; k < size; k++)
				MatrixResult[i][j] += MatrixA[i][k]* MatrixB[k][j];
		}
	}
}


//矩阵赋初值
void Matrix::MakeMatrix(long long int**  MatrixA, long long int size) { 
	for(long long int i = 0; i < size; i++) {
		for(long long int j = 0; j < size; j++) {
			MatrixA[i][j] = random(1,2);
		}
	}
}


void Matrix::StrassenMul(long long int**  MatrixA, long long int** MatrixB, long long int** MatrixResult, long long int size) {
	//矩阵阶数小于64，用传统矩阵乘法 
	if (size <= 64) {
		Matrix_Mul(MatrixA, MatrixB, MatrixResult, size);
		return; 
	}
	if(size == 1) {
		MatrixResult[0][0] = MatrixA[0][0]* MatrixB[0][0];
	} else {
		long long int** A11 = new long long int* [size/2];
		long long int** A12 = new long long int* [size/2];
		long long int** A21 = new long long int* [size/2];
		long long int** A22 = new long long int* [size/2];
		long long int** B11 = new long long int* [size/2];
		long long int** B12 = new long long int* [size/2];
		long long int** B21 = new long long int* [size/2];
		long long int** B22 = new long long int* [size/2];
		long long int** C11 = new long long int* [size/2];
		long long int** C12 = new long long int* [size/2];
		long long int** C21 = new long long int* [size/2];
		long long int** C22 = new long long int* [size/2];
		long long int** M1 = new long long int* [size/2];
		long long int** M2 = new long long int* [size/2];
		long long int** M3 = new long long int* [size/2];
		long long int** M4 = new long long int* [size/2];
		long long int** M5 = new long long int* [size/2];
		long long int** M6 = new long long int* [size/2];
		long long int** M7 = new long long int* [size/2];
		long long int** MatrixTemp1 = new long long int* [size/2];
		long long int** MatrixTemp2 = new long long int* [size/2];

		for(long long int i = 0; i < size/2; i++) {
			A11[i] = new long long int[size/2];
			A12[i] = new long long int[size/2];
			A21[i] = new long long int[size/2];
			A22[i] = new long long int[size/2];

			B11[i] = new long long int[size/2];
			B12[i] = new long long int[size/2];
			B21[i] = new long long int[size/2];
			B22[i] = new long long int[size/2];

			C11[i] = new long long int[size/2];
			C12[i] = new long long int[size/2];
			C21[i] = new long long int[size/2];
			C22[i] = new long long int[size/2];

			M1[i] = new long long int[size/2];
			M2[i] = new long long int[size/2];
			M3[i] = new long long int[size/2];
			M4[i] = new long long int[size/2];
			M5[i] = new long long int[size/2];
			M6[i] = new long long int[size/2];
			M7[i] = new long long int[size/2];

			MatrixTemp1[i] = new long long int[size/2];
			MatrixTemp2[i] = new long long int[size/2];
		}

		//赋值
		for(long long int i = 0; i < size/2; i++) {
			for(long long int j = 0; j < size/2; j++) {
				A11[i][j] = MatrixA[i][j];
				A12[i][j] = MatrixA[i][j+size/2];
				A21[i][j] = MatrixA[i+size/2][j];
				A22[i][j] = MatrixA[i+size/2][j+size/2];

				B11[i][j] = MatrixB[i][j];
				B12[i][j] = MatrixB[i][j+size/2];
				B21[i][j] = MatrixB[i+size/2][j];
				B22[i][j] = MatrixB[i+size/2][j+size/2];
			}
		}

		//M1
		Matrix_Sum(A11, A22, MatrixTemp1, size/2);
		Matrix_Sum(B11, B22, MatrixTemp2, size/2);
		StrassenMul(MatrixTemp1, MatrixTemp2, M1,size/2);

		//M2
		Matrix_Sum(A21, A22, MatrixTemp1, size/2);
		StrassenMul(MatrixTemp1, B11, M2, size/2);

		//M3
		Matrix_Sub(B12, B22, MatrixTemp1, size/2);
		StrassenMul(A11, MatrixTemp1, M3, size/2);


		//M4
		Matrix_Sub(B21, B11, MatrixTemp1, size/2);
		StrassenMul(A22, MatrixTemp1, M4, size/2);

		//M5
		Matrix_Sum(A11, A12, MatrixTemp1, size/2);
		StrassenMul(MatrixTemp1, B22, M5, size/2);

		//M6
		Matrix_Sub(A21, A11, MatrixTemp1, size/2);
		Matrix_Sum(B11, B12, MatrixTemp2, size/2);
		StrassenMul(MatrixTemp1, MatrixTemp2, M6, size/2);

		//M7
		Matrix_Sub(A12, A22, MatrixTemp1, size/2);
		Matrix_Sum(B21, B22, MatrixTemp2, size/2);
		StrassenMul(MatrixTemp1, MatrixTemp2, M7, size/2);

		//C11
		Matrix_Sum(M1, M4, C11, size/2);
		Matrix_Sub(C11, M5, C11, size/2);
		Matrix_Sum(C11, M7, C11, size/2);

		//C12
		Matrix_Sum(M3, M5, C12, size/2);

		//C21
		Matrix_Sum(M2, M4, C21, size/2);

		//C22
		Matrix_Sub(M1, M2, C22, size/2);
		Matrix_Sum(C22, M3, C22, size/2);
		Matrix_Sum(C22, M6, C22, size/2);

		//赋值
		for(long long int i = 0; i < size/2; i++) {
			for(long long int j = 0; j < size/2; j++) {
				MatrixResult[i][j] = C11[i][j];
				MatrixResult[i][j+size/2] = C12[i][j];
				MatrixResult[i+size/2][j] = C21[i][j];
				MatrixResult[i+size/2][j+size/2] = C22[i][j];
			}
		}

		//释放内存
		for(long long int i = 0; i < size/2; i++) {
			delete[] A11[i];
			delete[] A12[i];
			delete[] A21[i];
			delete[] A22[i];

			delete[] B11[i];
			delete[] B12[i];
			delete[] B21[i];
			delete[] B22[i];

			delete[] C11[i];
			delete[] C12[i];
			delete[] C21[i];
			delete[] C22[i];

			delete[] M1[i];
			delete[] M2[i];
			delete[] M3[i];
			delete[] M4[i];
			delete[] M5[i];
			delete[] M6[i];
			delete[] M7[i];

			delete[] MatrixTemp1[i];
			delete[] MatrixTemp2[i];
		}
		delete[] A11;
		delete[] A12;
		delete[] A21;
		delete[] A22;

		delete[] B11;
		delete[] B12;
		delete[] B21;
		delete[] B22;

		delete[] C11;
		delete[] C12;
		delete[] C21;
		delete[] C22;

		delete[] M1;
		delete[] M2;
		delete[] M3;
		delete[] M4;
		delete[] M5;
		delete[] M6;
		delete[] M7;

		delete[] MatrixTemp1;
		delete[] MatrixTemp2;
	}
}

//分块矩阵，分割矩阵A，将第x行第y列的元素作为新的矩阵元素，分割出size阶矩阵
Matrix Matrix::Matrix_Split(long long int** MatrixA, long long int x, long long int y, long long int size) {
	Matrix MatrixC(size);
	for (long long int i = 0; i < size; i++) {
		for (long long int j = 0; j < size; j++) {
			MatrixC.m[i][j] = MatrixA[i + x][j + y];
		}
	}
	return MatrixC;
}

//为阶数为偶数时的 Strassen算法
void Matrix::EvenNumberMul(long long int** matrixA , long long int** matrixB , long long int** MatrixResult , long long int size) {
	long long int n = size;
	long long int m = 1;
	//将矩阵拆为m*n,n为奇数，m为2的幂次
	//例如N=6;6=3*2^1
	while (!(n & 1)) {
		n >>= 1;//比特右移（>>）运算符,将 11100011 右移 3 比特，算术右移后成为 11111100，逻辑右移则为 00011100
		m <<= 1;
	}
	//将A拆分成n^2个m阶矩阵，在arr_A中存储,即用二维数组arr_A中的每一项[i][j]存储矩阵

	Matrix** arr_A = new Matrix* [n];
	for (long long int i = 0; i < n; i++) {
		arr_A[i] = new Matrix [n];//分配内存

		for (long long int j = 0; j < n; j++) {
			arr_A[i][j] = Matrix_Split(matrixA, i*m, j*m, m);
			//Matrix_Print(arr_A[i][j].m, m); //分割出来的矩阵
		}
	}
	//将B拆分成n^2个m阶矩阵，在arr_B中存储
	Matrix** arr_B = new Matrix*[n];
	for (long long int i = 0; i < n; i++) {
		arr_B[i] = new Matrix[n];//分配内存

		for (long long int j = 0; j < n; j++) {
			arr_B[i][j] = Matrix_Split(matrixB, i*m, j*m, m);
		}
	}
	//初始化数组arr_C ,存储arr_A *arr_B 的结果
	Matrix** arr_C = new Matrix* [n];
	for (long long int i = 0; i < n; i++) {
		arr_C[i] = new Matrix [n];//分配内存
		for (long long int j = 0; j < n; j++) {
			arr_C[i][j] = Matrix(m);
		}
	}

	//arr_C = arr_A * arr_B
	long long int** arr_temp = new long long int* [n];
	for(long long int t = 0; t < n; t++) arr_temp[t] = new long long int[n];
	for (long long int i = 0; i < n; i++) {
		for (long long int j = 0; j < n; j++) {
			for (long long int k = 0; k < n; k++) {
				StrassenMul(arr_A[i][k].m, arr_B[k][j].m, arr_temp, m);
				Matrix_Sum(arr_C[i][j].m, arr_temp, arr_C[i][j].m, m);
			}
		}
	}
	//声明矩阵C(m*n=N)，将各矩阵合并为一个矩阵
	for (long long int i = 0; i < n; i++) { //合并
		for (long long int j = 0; j < n; j++) {
			for (long long int x = 0; x < m; x++) {
				for (long long int y = 0; y < m; y++) {
					MatrixResult[x + i*m][y + j*m] = arr_C[i][j].m[x][y];
				}
			}
		}
	}
//	Matrix_Print(MatrixResult, m*n);
	//释放申请的内存
	for(long long int i = 0; i < n; i++) {
		delete[] arr_A[i];
		delete[] arr_B[i];
		delete[] arr_C[i];
	}
	delete[] arr_A;
	delete[] arr_B;
	delete[] arr_C;
}

//为阶数为奇数时的 Strassen算法
void Matrix::OddNumberMul(long long int** matrixA,long long int** matrixB,long long int** MatrixResult,long long int size) {
	//扩容矩阵，第一行，第一列增加，A[0,0] 位置的值为1
	size=size+1;
	long long int** newA=new long long int* [size];
	long long int** newB=new long long int* [size];
	long long int** newC=new long long int* [size];

	for(long long int i=0; i<size; i++) {
		newA[i]=new long long int[size];
		newB[i]=new long long int[size];
		newC[i]=new long long int[size];
	}

	for(long long int i=0; i<size; i++) {
		for(long long int j=0; j<size; j++) {
			if(i==0||j==0) {
				newA[i][j]=0;
				newB[i][j]=0;
				continue;
			}
			newA[i][j]=matrixA[i-1][j-1];
			newB[i][j]=matrixB[i-1][j-1];
		}
	}

	newA[0][0]=1;
	newB[0][0]=1;
	if(IsPower(size) == 1) {
		StrassenMul(newA, newB, newC, size);
	} else {
		EvenNumberMul(newA, newB, newC, size);
	}
	for(long long int i=0; i<size-1; i++) {
		for (long long int j = 0; j < size-1; j++) {
			MatrixResult[i][j] = newC[i+1][j+1];
		}
	}
	//释放申请的内存
	for(long long int i = 0; i < size; i++) {
		delete[] newA[i];
		delete[] newB[i];
		delete[] newC[i];
	}
	delete[] newA;
	delete[] newB;
	delete[] newC;
}


//综合所有情况
void Matrix::Matrix_Mul_AllType(long long int** matrixA,long long int** matrixB,long long int** matrixC,long long int size) {
	if(size % 2 == 1) {
		cout<<"阶数为奇数：\n";
		OddNumberMul(matrixA, matrixB, matrixC, size);
	} else {
		if(IsPower(size) == 1) {
			cout<<"阶数为 2 的 K 次方：\n";
			StrassenMul(matrixA, matrixB, matrixC, size);
		} else {
			cout<<"阶数为非 2 的幂次，且为偶数：\n";
			EvenNumberMul(matrixA, matrixB, matrixC, size);
		}
	}
}


long long int Matrix::GetMatrixSum(long long int** Matrix, long long int size) {
	long long int sum = 0;
	for(long long int i = 0; i < size; i++) {
		for(long long int j = 0; j < size; j++) {
			sum += Matrix[i][j];
		}
	}
	return sum;
}

int main() {
	//srand(time(0));
	Matrix matrix;
	long long int startTime_normal, endTime_normal;
	long long int startTime_strasse, endTime_strassen;
	long long int N;
	cout<<"输入矩阵阶数：";
	cin>>N;
	cout<<endl; 

	long long int** Matrix1 = new long long int* [N];
	long long int** Matrix2 = new long long int* [N];
	long long int** Matrix3 = new long long int* [N];
	for(long long int i=0; i<N; i++) {
		Matrix1[i] = new long long int[N];
		Matrix2[i] = new long long int[N];
		Matrix3[i] = new long long int[N];
	}

	matrix.MakeMatrix(Matrix1, N);
	matrix.MakeMatrix(Matrix2, N);
//	cout<<"随机生成的矩阵1:\n";
//	matrix.Matrix_Print(Matrix1, N);
//	cout<<"随机生成的矩阵2:\n";
//	matrix.Matrix_Print(Matrix2, N);

	cout<<"传统算法开始时间："<<(startTime_normal = clock())<<endl;
	matrix.Matrix_Mul(Matrix1, Matrix2, Matrix3, N);
	cout<<"传统算法结束时间："<<(endTime_normal = clock())<<endl;
	cout<<"总时间:"<<endTime_normal-startTime_normal<<endl;
	cout<<"sum = "<<matrix.GetMatrixSum(Matrix3, N)<<endl;
//	cout<<"传统算法结果:\n";
//	matrix.Matrix_Print(Matrix3, N);
	cout<<"\n--------------------------------\n"; 
	cout<<"Strassen算法开始时间："<<(startTime_strasse= clock())<<endl;
	matrix.Matrix_Mul_AllType(Matrix1, Matrix2, Matrix3, N);
	cout<<"Strassen算法结束时间："<<(endTime_strassen = clock())<<endl;
	cout<<"总时间:"<<endTime_strassen-startTime_strasse<<endl;
	cout<<"sum = "<<matrix.GetMatrixSum(Matrix3, N)<<endl;
//	cout<<"Strassen算法结果:\n";
//	matrix.Matrix_Print(Matrix3, N);
	return 0;
}

