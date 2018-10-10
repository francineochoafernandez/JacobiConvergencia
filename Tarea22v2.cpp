#include <iostream>
#include <string>
#include <cstring>
#include <cmath>

using namespace std;

int i,j,k;

typedef struct Jacobi
{
  float a[20][20],br[20][20], x[20],b[20],trInf[20][20],trSup[20][20],DiaInv[20][20],MatT[20][20],VectC[20],MatAux1[20][20],VectAux1[20],VectAux2[20],VectAux3[20];
  int n, maxit, indi;
  float err;

  void RellenaMat(void)
  {
    printf("\n\nDar la cantidad de incógnitas de la ecuación: ");
    scanf("%d",&n);
    for(i=0;i<n;i++)
    {
        printf("\n\nDa todos los coeficientes de la ecuacion %d: \n",i+1);
        for(j=0;j<n;j++)
          scanf("%f",&a[i][j]);
    }

    printf("\n\nDa los coeficientes del vector b: \n");
    for(j=0;j<n;j++)
      scanf("%f",&b[j]);

    printf("\n\nDar el vector x con la estimación de los resultados: \n");
    for(j=0;j<n;j++)
      scanf("%f",&x[j]);

    printf("\n\nDar el error relativo y el numero de iteraciones:  \n");
    scanf("%f%d",&err,&maxit);

  }

  void ImprimeOp(void)
  {
    //Imprimiendo matrices
    for(i=0;i<n;i++)
    {
      cout << "|";
      for(j=0;j<n;j++)
      {
        printf("\t%f \t",a[i][j]);
      }
      cout << "|";
      if(i==n/2)
      {
        printf("   *   ");
      }
      else
      {
        printf("       ");
      }

      printf("| %f |\t",x[i]);

      if(i==n/2)
      {
        printf("   =   ");
      }
      else
      {
        printf("       ");
      }

      printf("| %f | \t",b[i]);
      printf("\n");
    }
    printf("\n");
  }

  void CopiaMatriz(float mat1[20][20],float mat2[20][20])
  {
    for(i=0;i<n;i++)
    {
      for(j=0;j<n;j++)
      {
        mat2[i][j]=mat1[i][j];
      }
    }
  }

  void AcomodandoMatriz()
  {
    float renc, renc2,aux;
    int y,r=0;

    //Combirtiendola en diagonal Dominante
    for(i=0;i<n;i++)
    {
      renc=0;
      for(j=0;j<n;j++)
      {
        if (i!=j)
        {
          renc= renc + abs(br[i][j]);
        }
      }

      y=0;
      r=i;
      if(renc>abs(br[i][i]))//El renglon no cumple así que se cambia orden
      {
        for(int l=0;l<n;l++)//Este for solo es para que se repitan las lineas la n cantidad de veces.
        {//Aquí se hacen todas las posibles combinaciones para ver cual digito es el que iria en la diagonal
          renc2=0;
          aux=0;
          for(int k=0;k<n;k++)
          {
            if(k!=y)
              renc2= renc2 + abs(br[i][k]);
          }

          if(abs(br[i][y])>renc2)
          {
            for(i=0;i<n;i++)
            {
  	           aux=br[i][r];
  	           br[i][r]=br[i][y];
  	           br[i][y]=aux;
            }
          }
          y++;

        }

      }

    }

  }

  void ComparaDiago(int op)
  {
    //Comparando diagonal
    float ren;
    int c=0;

    for(i=0;i<n;i++)
    {
      ren=0;
      for(j=0;j<n;j++)
      {
        if (i!=j)
        {
          ren= ren + abs(br[i][j]);
        }
      }

      if(ren>abs(br[i][i]))
      {
        if(op==1)
        {
          printf("\nLa matriz NO es dominante, no va a convergir.\n");
          printf("\nSe tratara de mover la matriz.\n");
          AcomodandoMatriz();
        }
        else
        {
          printf("\nLa matriz NO se pudo acomodar para volverla dominante.\n");
          indi=3;
        }
        i=n;
        j=n;
        c++;
      }

    }
    if(c==0)
    {
      printf("\nLa matriz SI es dominante, sí va a convergir.\n");
      indi=1;
      if(op==2)
      {
        CopiaMatriz( br, a);
        ImprimeOp();
      }

    }
  }



  void MatricesTriangulares()
  { //Saca la matriz triangular inferior a partir de a
    for(i=0;i<n;i++)
    {
      for (j=0;j<n;j++)
      {
        if (j>i)
        {
          trSup[i][j]=a[i][j];
        }
        else
        {
          trSup[i][j]=0;
        }


      }
    }

    for(i=0;i<n;i++)
    {
      for (j=0;j<n;j++)
      {
        if (j<i)
        {
          trInf[i][j]=a[i][j];
        }
        else
        {
          trInf[i][j]=0;
        }
      }
    }

  }

  void ImprimeMat(float matriz[20][20])
  {
    for(i=0;i<n;i++)
    {
      cout << "|";
      for(j=0;j<n;j++)
      {
        printf("\t %f \t",matriz[i][j]);
      }
      cout << "|";
      printf("\n");
    }
    printf("\n");
  }

  void Mult(float matriz1[20][20],float matriz2[20][20], float matrizR[20][20])
  {
    float sum=0.0;
    for (i=0;i<n;i++)
    {
      for (j=0;j<n;j++)
      {
        for (int k=0;k<n;k++)
        {
          sum = sum + matriz1[i][k]*matriz2[k][j];
        }
        matrizR[i][j] = sum;
        sum = 0.0;
      }
    }
  }

  void MultVec(float matriz1[20][20],float vect[20], float vectR[20])
  {
    float sum=0.0;
    for (i=0;i<n;i++)
    {
      for (j=0;j<n;j++)
      {
        sum = sum + matriz1[i][j]*vect[j];
        vectR[i] = sum;
      }
      sum = 0.0;
    }
  }

  void Suma(float matriz1[20][20],float matriz2[20][20], float matrizR[20][20])
  {
    //Sumando matrices
    for(int i=0;i<n;i++)
    {
      for(int j=0;j<n;j++)
      {
        matrizR[i][j]=matriz1[i][j]+matriz2[i][j];
      }
    }
  }

  void RestVect(float vec1[20],float vec2[20], float vecR[20])
  {
    for(int i=0;i<n;i++)
    {
      vecR[i]=vec1[i]-vec2[i];
    }
  }

  void SumaVect(float vec1[20],float vec2[20], float vecR[20])
  {
    //Sumando matrices
    for(int i=0;i<n;i++)
    {
      vecR[i]=vec1[i]+vec2[i];
    }
  }

  void DiagonalInversa()
  {
    for(i=0;i<n;i++)
    {
      for(j=0;j<n;j++)
      {
        if (i==j)
        {
          DiaInv[i][j]= 1/(a[i][j]);
        }
        else
        {
          DiaInv[i][j]=0;
        }
      }
    }
  }

  void MatrizTyC()//T=-D^(-1) * (trInf+trSup)
  {
    MatricesTriangulares();//Se sacan matrices triangulares
    DiagonalInversa();//Se hace la matriz de la diagonal inversa
    Suma(trInf,trSup,MatAux1);
    Mult(DiaInv,MatAux1,MatT);
    //Multiplicando por -1
    for(i=0;i<n;i++)
    {
      for(j=0;j<n;j++)
      {
        MatT[i][j]=(-1)*MatT[i][j];
      }
    }
    MultVec(DiaInv,b,VectC);
  }

  float Criterio()
  {
    float mag=0.0;
    MultVec(a,x,VectAux2);
    RestVect(VectAux2,b,VectAux3);
    for(i=0;i<n;i++)
    {
      mag= mag+ VectAux3[i]*VectAux3[i];
    }
    return sqrt(mag);
  }
  void Iteraciones ()
  {
    int c=0;
    for(i=0;i<maxit;i++)
    {
      MultVec(MatT,x,VectAux1);
      SumaVect(VectAux1,VectC,x);

      if(Criterio()<=err)
      {
        printf("\nEl vector SOLUCION es:\n");
        for(i=0;i<n;i++)
          printf("| %f |\n",x[i]);
        printf("\n La solución se obtuvo en %i iteraciones.\n",c);
        i=maxit;
      }

      c++;
    }

  }


}Jac;


int main(void)
{
  Jac ecs;
  ecs.RellenaMat();//Se piden valores de matriz y demas
  ecs.ImprimeOp();//Se imprimen los valores dados
  ecs.CopiaMatriz( ecs.a, ecs.br);
  ecs.ComparaDiago(1);
  if(ecs.indi!=1)
  {
    ecs.ComparaDiago(2);
  }

  ecs.MatrizTyC();
  /*
  printf("\nMatriz triangular Inferior:\n");
  ecs.ImprimeMat(ecs.trInf);
  printf("\nMatriz triangular Superior:\n");
  ecs.ImprimeMat(ecs.trSup);
  printf("\nSuma de las matrices triangulares Inferior y Superior:\n");
  ecs.ImprimeMat(ecs.MatAux1);
  printf("\nMatriz con diagonal inversa:\n");
  ecs.ImprimeMat(ecs.DiaInv);
  */
  if(ecs.indi!=3)
  {
    ecs.Iteraciones();
  }

  return 0;
}
//https://www.wikiwand.com/es/M%C3%A9todo_de_Jacobi
