// FFT.cpp : 定义控制台应用程序的入口点。
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
//定义最大序列长度
const int N = 2048;
int n;
const float PI = 3.1416;



/////快速傅里叶变换部分
inline void swap(float &a, float &b)
{
	float t;
	t = a;
	a = b;
	b = t;
}

void bitrp(float xreal[], float ximag[], int n)
{
	// 位反转置换
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
	// 快速傅立叶变换，将复数 x 变换后仍保存在 x 中，xreal, ximag 分别是 x 的实部和虚部
	float wreal[N / 2], wimag[N / 2], treal, timag, ureal, uimag, arg;
	int m, k, j, t, gIn, hIn;

	bitrp(xreal, ximag, n);

	// 计算 前 n / 2 个 n 次方根的共轭复数 W'j = wreal [j] + i * wimag [j] , j = 0, 1, ... , n / 2 - 1
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
				t = n * j / m;    // 旋转因子 w 的实部在 wreal [] 中的下标为 t
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
	// 快速傅立叶反变换
	float wreal[N / 2], wimag[N / 2], treal, timag, ureal, uimag, arg;
	int m, k, j, t, gIn, hIn;

	bitrp(xreal, ximag, n);

	// 计算 1 的前 n / 2 个 n 次方根 Wj = wreal [j] + i * wimag [j] , j = 0, 1, ... , n / 2 - 1
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
				t = n * j / m;    // 旋转因子 w 的实部在 wreal [] 中的下标为 t
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
//////////////////////////快速傅里叶变换结束
///////////////////////调用OpenGL显示
void init()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void reshape(int w, int h)
{
	glClearColor(1.0, 1.0, 1.0, 1.0); //设置背景为白色
	glViewport(0, 0, w, h);
	glLoadIdentity();
	gluPerspective(45.0, (GLfloat)w / (GLfloat)h, -1.0, 20.0);
	glMatrixMode(GL_PROJECTION);



}

void myDisplay(void)
{

	glClearColor(1.0, 1.0, 1.0, 1.0); //设置背景为白色
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
	//	左上，三角函数的原始图像

	glPushMatrix();
	glTranslated(-0.5, 0.5, 0);
	glScaled(0.5, 0.5, 0.5);
	//绘制坐标系

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

	//绘制点，并对原始的点赋值
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
	cout << "三角形函数对应的值" << endl;
	cout << "下标     实部    虚部     " << endl;

	for (int i = 0;i < n;i++)
	{

		cout << i << "       " << xreal[i] << "         " << ximag[i] << endl;

	}
	FFT(xreal, ximag, n);
	cout << "傅里叶变换对应的值" << endl;
	cout << "下标        实部         虚部     " << endl;
	for (int i = 0;i < n;i++)
	{
		cout << i << "         " << xreal[i] << "           " << ximag[i] << endl;
	}

	//右上，三角函数的傅里叶变换
	glPushMatrix();
	glTranslated(0.5, 0.5, 0);
	glScaled(0.5, 0.5, 0.5);
	//绘制坐标系

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
	//黑色表示实部
	glColor3f(0, 0, 0);
	glBegin(GL_LINES);
	for (int i = 1;i <= n;i++)
	{
		glVertex3f(double(i) / n, 0, 0);
		glVertex3f(double(i) / n, xreal[i]/n, 0);
	}
	glEnd();
	//蓝色表示虚部
	glColor3f(0, 0,1);
	glBegin(GL_LINES);
	for (int i = 1;i < n;i++)
	{
		glVertex3f(double(i) /n,0, 0);
		glVertex3f(double(i) / n, ximag[i] / (5*n), 0);
	}
	glEnd();

	IFFT(xreal, ximag, n);
	//红色显示傅里叶反变换的结果
	glColor3d(1, 0, 0);
	glBegin(GL_POINTS);
	for (int i = 1;i < n;i++)
	{
		glVertex3f(double(i) / n, xreal[i]/n, 0);
	}
	glEnd();




	glPopMatrix();
	//左下，正弦函数的原始图像
	glPushMatrix();
	glTranslated(-0.5, -0.5, 0);
	glScaled(0.5, 0.5, 0.5);
	//绘制坐标系
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
	cout << "正弦函数对应的值" << endl;
	cout << "下标     实部    虚部     " << endl;

	for (int i = 0;i < n;i++)
	{

		cout << i << "       " << xreal[i] << "         " << ximag[i] << endl;

	}
	//进行傅里叶变换
	FFT(xreal, ximag, n);
	cout << "傅里叶变换对应的值" << endl;
	cout << "下标        实部         虚部     " << endl;
	for (int i = 0;i < n;i++)
	{
		cout << i << "         " << xreal[i] << "           " << ximag[i] << endl;
	}

	//右下，正弦函数的快速傅里叶变换结果
	glPushMatrix();
	glTranslated(0.5,- 0.5, 0);
	glScaled(0.5, 0.5, 0.5);
	//绘制坐标系
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

	//黑色表示实部
	glColor3f(0, 0, 0);
	glBegin(GL_LINES);
	for (int i = 1;i <= n;i++)
	{
		glVertex3f(double(i) / n, 0, 0);
		glVertex3f(double(i) / n, xreal[i] / n, 0);
	}
	glEnd();

	//蓝色表示虚部
	glColor3f(0, 0, 1);
	glBegin(GL_LINES);
	for (int i = 1;i < n;i++)
	{
		glVertex3f(double(i) / n, 0, 0);
		glVertex3f(double(i) / n, ximag[i]/n , 0);
	}
	glEnd();

	IFFT(xreal, ximag, n);
	//红色显示傅里叶反变换的结果
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

	cout << "请输入长度m（序列长度N = 2 ^ m, 为了显示效果请控制m <= 10!" << endl;
	cin >> m;
	n = pow(2, m);
	glutCreateWindow("FFT变换");
	glutDisplayFunc(myDisplay);
	glutReshapeFunc(reshape);

	glutMainLoop();
	return 0;
}
