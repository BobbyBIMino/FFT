// FFT.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
#include<GL/glut.h>
#include <math.h>
#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<iostream>

using namespace std;
int m;
//����������г���
const int N = 2048;
int n;
const float PI = 3.1416;



/////���ٸ���Ҷ�任����
inline void swap(float &a, float &b)
{
	float t;
	t = a;
	a = b;
	b = t;
}

void bitrp(float xreal[], float ximag[], int n)
{
	// λ��ת�û�
	int i, j, a, b, p;

	for (i = 1, p = 0; i < n; i *= 2)
	{
		p++;
	}
	for (i = 0; i < n; i++)
	{
		a = i;
		b = 0;
		for (j = 0; j < p; j++)
		{
			b = (b << 1) + (a & 1);    // b = b * 2 + a % 2;
			a >>= 1;        // a = a / 2;
		}
		if (b > i)
		{
			swap(xreal[i], xreal[b]);
			swap(ximag[i], ximag[b]);
		}
	}
}

void FFT(float xreal[], float ximag[], int n)
{
	// ���ٸ���Ҷ�任�������� x �任���Ա����� x �У�xreal, ximag �ֱ��� x ��ʵ�����鲿
	float wreal[N / 2], wimag[N / 2], treal, timag, ureal, uimag, arg;
	int m, k, j, t, gIn, hIn;

	bitrp(xreal, ximag, n);

	// ���� ǰ n / 2 �� n �η����Ĺ���� W'j = wreal [j] + i * wimag [j] , j = 0, 1, ... , n / 2 - 1
	arg = -2 * PI / n;
	treal = cos(arg);
	timag = sin(arg);
	wreal[0] = 1.0;
	wimag[0] = 0.0;
	for (j = 1; j < n / 2; j++)
	{
		wreal[j] = wreal[j - 1] * treal - wimag[j - 1] * timag;
		wimag[j] = wreal[j - 1] * timag + wimag[j - 1] * treal;
	}

	for (m = 2; m <= n; m *= 2)
	{
		for (k = 0; k < n; k += m)
		{
			for (j = 0; j < m / 2; j++)
			{
				gIn = k + j;
				hIn = gIn + m / 2;
				t = n * j / m;    // ��ת���� w ��ʵ���� wreal [] �е��±�Ϊ t
				treal = wreal[t] * xreal[hIn] - wimag[t] * ximag[hIn];
				timag = wreal[t] * ximag[hIn] + wimag[t] * xreal[hIn];
				ureal = xreal[gIn];
				uimag = ximag[gIn];
				xreal[gIn] = ureal + treal;
				ximag[gIn] = uimag + timag;
				xreal[hIn] = ureal - treal;
				ximag[hIn] = uimag - timag;
			}
		}
	}
}

void  IFFT(float xreal[], float ximag[], int n)
{
	// ���ٸ���Ҷ���任
	float wreal[N / 2], wimag[N / 2], treal, timag, ureal, uimag, arg;
	int m, k, j, t, gIn, hIn;

	bitrp(xreal, ximag, n);

	// ���� 1 ��ǰ n / 2 �� n �η��� Wj = wreal [j] + i * wimag [j] , j = 0, 1, ... , n / 2 - 1
	arg = 2 * PI / n;
	treal = cos(arg);
	timag = sin(arg);
	wreal[0] = 1.0;
	wimag[0] = 0.0;
	for (j = 1; j < n / 2; j++)
	{
		wreal[j] = wreal[j - 1] * treal - wimag[j - 1] * timag;
		wimag[j] = wreal[j - 1] * timag + wimag[j - 1] * treal;
	}

	for (m = 2; m <= n; m *= 2)
	{
		for (k = 0; k < n; k += m)
		{
			for (j = 0; j < m / 2; j++)
			{
				gIn = k + j;
				hIn = gIn + m / 2;
				t = n * j / m;    // ��ת���� w ��ʵ���� wreal [] �е��±�Ϊ t
				treal = wreal[t] * xreal[hIn] - wimag[t] * ximag[hIn];
				timag = wreal[t] * ximag[hIn] + wimag[t] * xreal[hIn];
				ureal = xreal[gIn];
				uimag = ximag[gIn];
				xreal[gIn] = ureal + treal;
				ximag[gIn] = uimag + timag;
				xreal[hIn] = ureal - treal;
				ximag[hIn] = uimag - timag;
			}
		}
	}

	for (j = 0; j < n; j++)
	{
		xreal[j] /= n;
		ximag[j] /= n;
	}
}
//////////////////////////���ٸ���Ҷ�任����
///////////////////////����OpenGL��ʾ
void init()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void reshape(int w, int h)
{
	glClearColor(1.0, 1.0, 1.0, 1.0); //���ñ���Ϊ��ɫ
	glViewport(0, 0, w, h);
	glLoadIdentity();
	gluPerspective(45.0, (GLfloat)w / (GLfloat)h, -1.0, 20.0);
	glMatrixMode(GL_PROJECTION);



}

void myDisplay(void)
{

	glClearColor(1.0, 1.0, 1.0, 1.0); //���ñ���Ϊ��ɫ
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(0.0, 0.0, 0.0);
	float xreal[N] = {}, ximag[N] = {};
	init();
	glLineWidth(5);
	glBegin(GL_LINES);
	glVertex3f(-1.0f, 0.0f, 0.0f);
	glVertex3f(1.0, 0.0f, 0.0f);
	glVertex3f(0.0, -1.0f, 0.0f);
	glVertex3f(0.0, 1.0f, 0.0f);
	glEnd();
	glLineWidth(2);
	glPointSize(20.0 / m);
	//	���ϣ����Ǻ�����ԭʼͼ��

	glPushMatrix();
	glTranslated(-0.5, 0.5, 0);
	glScaled(0.5, 0.5, 0.5);
	//��������ϵ

	glBegin(GL_LINES);
	glVertex3f(-1.0f, 0.0f, 0.0f);
	glVertex3f(1.0, 0.0f, 0.0f);
	glVertex3f(0.0, -1.0f, 0.0f);
	glVertex3f(0.0, 1.0f, 0.0f);
	for (double i = -10;i <= 10;i++)
	{
		glVertex3f(i / 10, 0, 0.0f);
		glVertex3f(i / 10, 0.03f, 0.0f);
	}
	glEnd();

	//���Ƶ㣬����ԭʼ�ĵ㸳ֵ
	glBegin(GL_POINTS);
	for (double i = 0;i < n / 4;i++)
	{
		xreal[int(i)] = i;
		ximag[int(i)] = 0;
		glVertex3f(i / n, i / n, 0.0f);
	}
	for (double i = n / 4;i < n / 4 * 3;i++)
	{
		xreal[int(i)] = (-i / n + 0.5)*n;
		ximag[int(i)] = 0;
		glVertex3f(i / n, -i / n + 0.5, 0.0f);
	}

	for (double i = n / 4 * 3;i < n;i++)
	{
		xreal[int(i)] = (i / n - 1)*n;
		ximag[int(i)] = 0;
		glVertex3f(i / n, i / n - 1, 0.0f);
	}
	glEnd();
	glPopMatrix();
	cout << "�����κ�����Ӧ��ֵ" << endl;
	cout << "�±�     ʵ��    �鲿     " << endl;

	for (int i = 0;i < n;i++)
	{

		cout << i << "       " << xreal[i] << "         " << ximag[i] << endl;

	}
	FFT(xreal, ximag, n);
	cout << "����Ҷ�任��Ӧ��ֵ" << endl;
	cout << "�±�        ʵ��         �鲿     " << endl;
	for (int i = 0;i < n;i++)
	{
		cout << i << "         " << xreal[i] << "           " << ximag[i] << endl;
	}

	//���ϣ����Ǻ����ĸ���Ҷ�任
	glPushMatrix();
	glTranslated(0.5, 0.5, 0);
	glScaled(0.5, 0.5, 0.5);
	//��������ϵ

	glBegin(GL_LINES);
	glVertex3f(-1.0f, 0.0f, 0.0f);
	glVertex3f(1.0, 0.0f, 0.0f);
	glVertex3f(0.0, -1.0f, 0.0f);
	glVertex3f(0.0, 1.0f, 0.0f);
	for (double i = -10;i <= 10;i++)
	{
		glVertex3f(i / 10, 0, 0.0f);
		glVertex3f(i / 10, 0.03f, 0.0f);
	}
	glEnd();
	//��ɫ��ʾʵ��
	glColor3f(0, 0, 0);
	glBegin(GL_LINES);
	for (int i = 1;i <= n;i++)
	{
		glVertex3f(double(i) / n, 0, 0);
		glVertex3f(double(i) / n, xreal[i]/n, 0);
	}
	glEnd();
	//��ɫ��ʾ�鲿
	glColor3f(0, 0,1);
	glBegin(GL_LINES);
	for (int i = 1;i < n;i++)
	{
		glVertex3f(double(i) /n,0, 0);
		glVertex3f(double(i) / n, ximag[i] / (5*n), 0);
	}
	glEnd();

	IFFT(xreal, ximag, n);
	//��ɫ��ʾ����Ҷ���任�Ľ��
	glColor3d(1, 0, 0);
	glBegin(GL_POINTS);
	for (int i = 1;i < n;i++)
	{
		glVertex3f(double(i) / n, xreal[i]/n, 0);
	}
	glEnd();




	glPopMatrix();
	//���£����Һ�����ԭʼͼ��
	glPushMatrix();
	glTranslated(-0.5, -0.5, 0);
	glScaled(0.5, 0.5, 0.5);
	//��������ϵ
	glColor3d(0, 0, 0);

	glBegin(GL_LINES);
	glVertex3f(-1.0f, 0.0f, 0.0f);
	glVertex3f(1.0, 0.0f, 0.0f);
	glVertex3f(0.0, -1.0f, 0.0f);
	glVertex3f(0.0, 1.0f, 0.0f);
	for (double i = -10;i <= 10;i++)
	{
		glVertex3f(i / 10, 0, 0.0f);
		glVertex3f(i / 10, 0.03f, 0.0f);
	}
	glEnd();
	glBegin(GL_POINTS);
	for (double i = 0;i <= n;i++)
	{
		glVertex3f(i / n, cos(i / n * 2 * PI), 0.0f);
		xreal[int(i)] = cos(i / n * 2 * PI);
		ximag[int(i)] = 0;
	}
	glEnd();

	glPopMatrix();
	cout << "���Һ�����Ӧ��ֵ" << endl;
	cout << "�±�     ʵ��    �鲿     " << endl;

	for (int i = 0;i < n;i++)
	{

		cout << i << "       " << xreal[i] << "         " << ximag[i] << endl;

	}
	//���и���Ҷ�任
	FFT(xreal, ximag, n);
	cout << "����Ҷ�任��Ӧ��ֵ" << endl;
	cout << "�±�        ʵ��         �鲿     " << endl;
	for (int i = 0;i < n;i++)
	{
		cout << i << "         " << xreal[i] << "           " << ximag[i] << endl;
	}

	//���£����Һ����Ŀ��ٸ���Ҷ�任���
	glPushMatrix();
	glTranslated(0.5,- 0.5, 0);
	glScaled(0.5, 0.5, 0.5);
	//��������ϵ
	glBegin(GL_LINES);
	glVertex3f(-1.0f, 0.0f, 0.0f);
	glVertex3f(1.0, 0.0f, 0.0f);
	glVertex3f(0.0, -1.0f, 0.0f);
	glVertex3f(0.0, 1.0f, 0.0f);
	for (double i = -10;i <= 10;i++)
	{
		glVertex3f(i / 10, 0, 0.0f);
		glVertex3f(i / 10, 0.03f, 0.0f);
	}
	glEnd();

	//��ɫ��ʾʵ��
	glColor3f(0, 0, 0);
	glBegin(GL_LINES);
	for (int i = 1;i <= n;i++)
	{
		glVertex3f(double(i) / n, 0, 0);
		glVertex3f(double(i) / n, xreal[i] / n, 0);
	}
	glEnd();

	//��ɫ��ʾ�鲿
	glColor3f(0, 0, 1);
	glBegin(GL_LINES);
	for (int i = 1;i < n;i++)
	{
		glVertex3f(double(i) / n, 0, 0);
		glVertex3f(double(i) / n, ximag[i]/n , 0);
	}
	glEnd();

	IFFT(xreal, ximag, n);
	//��ɫ��ʾ����Ҷ���任�Ľ��
	glColor3d(1, 0, 0);
	glBegin(GL_POINTS);
	for (int i = 1;i <= n;i++)
	{
		glVertex3f(double(i) / n, xreal[i] , 0);
	}
	glEnd();
	glPopMatrix();

	glFlush();
}




int main(int argc, char *argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
	glutInitWindowPosition(100, 100);
	gluOrtho2D(0, 1200, 0, 1200);
	glutInitWindowSize(1000, 1000);

	cout << "�����볤��m�����г���N = 2 ^ m, Ϊ����ʾЧ�������m <= 10!" << endl;
	cin >> m;
	n = pow(2, m);
	glutCreateWindow("FFT�任");
	glutDisplayFunc(myDisplay);
	glutReshapeFunc(reshape);

	glutMainLoop();
	return 0;
}
