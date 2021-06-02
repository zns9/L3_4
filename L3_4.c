//��ᐧ��ɂ�鑬�x����v���O���� L3_4.c�i���ӁF�v���O������̃X�y�[�X�͑S�āu���p�v

#include <dos.h>
#include <math.h>
#include <stdio.h>
#include <conio.h>
#include "shfunc.h"

int main(void)
{
	FILE *fp1;                                                              // �t�@�C���|�C���^�̐ݒ�
	float x1, x2, x3, x4, x5, x6, x7, x8, x9;                    // �g�p����ϐ��̒�`
	float c1, c2, c3, c4, c5;
	float Ra, La, K, J, D, Va, TL, TLs, TLe;
	float ts, te, dt, taue, tauj;
	float t, Tload;
	float A, Tau1, Vmin, Vmax, TG, Wcmd;

	x4 = x8 = x9 = c5 = 0.0;                                       // �g�p����ϐ��̏�����

	Ra = 1.0;                                                               // �d�@�q��R
	La = 0.008;                                                           // �d�@�q�C���_�N�^���X

	K = 0.28;                                                               // ���[�^�萔�iKE = kT = K�j
	J = 0.005;                                                              // ���ׂ��܂ފ������[�����g
	D = 0.0002;                                                           // �S�������W��

	Wcmd = 5.0;                                                         // ���x�w�߁i1500/min�j
	Vmin = -80.0;                                                        // �ŏ����[�^�쓮�d��
	Vmax = 80.0;                                                         // �ő僂�[�^�쓮�d��
	TL = 1.4;                                                               // ���׃g���N
	ts = 0.0;                                                                // �V�~�����[�V�����J�n����
	te = 2.0;                                                                // �V�~�����[�V�����I������
	TLs = 0.4;                                                              // ���׈���J�n����
	TLe = 0.7;                                                              // ���׈���I������

	dt = 0.0010;                                                           // �v�Z���ԊԊu
	taue = La / Ra;
	tauj = J / D;
	Tau1 = 0.01;                                                           // ���[�p�X�t�B���^�̐ݒ�
	A = 76.25;                                                                // ���x�A���v�̃Q�C����ݒ�
	TG = 10.0 / 314.159;                                               // ���x���o��̐ݒ�
                                                                                             // 3000/min �̎�, 10V ���o��

    fp1 = fopen("c:\\Users\\yoshi\\data4.76.25.csv", "w");//作成するファイルのオープン
	fprintf(fp1, "%f, %f, %f\n", ts, 0.0, 0.0);                    // �����f�[�^�̃t�@�C���ւ̏�������

	for(t = dt; t < te; t += dt)                                        // dt�`te�b�܂ł̌J��Ԃ�
	{
		if(t > TLs && t < TLe) Tload = TL; else Tload = 0.0;            
                                                                                              // TLs < t < TLe �̎��ATload = TL
                                                                                              // ����ȊO�́ATload = 0.0
		c1 = Wcmd;                                              // ���x�w�߂̐ݒ�
		c2 = c1 - c5;
		c3 = c2 * A;
		x1 = limit(c3, Vmax, Vmin);                       // �d�@�q�d���̐���
		x2 = x1 - x9;
		x3 = x2 / Ra;
		x4 = rk4(x3, x4, taue, dt);                         // �d�@�q�d���̌v�Z
		x5 = x4 * K;                                              // ���[�^�����g���N�̌v�Z
		x6 = x5 - Tload;                                         // ���׃g���N�����
		x7 = x6 / D;
		x8 = rk4(x7, x8, tauj, dt);                          // ��]�p���x�̌v�Z
		x9 = x8 * K;                                              // �d���̌v�Z
		c4 = x8 * TG;                                            // ��]���x�̌��o
                              c5 = rk4(c4, c5, Tau1, dt);                         // �m�C�Y�����p���[�p�X�t�B���^

		fprintf(fp1, "%f, %f, %f\n", t, x8 * 60.0 / 6.28318, x4);           
                                                                   // ����(sec), �d��(A), ��]�p���x(min^-1), ���׃g���N
                                                                   // �f�[�^�̃t�@�C���ւ̏�������
	}
	fclose(fp1);                                                            // �t�@�C���̃N���[�Y
	return(0);
}