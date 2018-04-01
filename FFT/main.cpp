//�˳�����ÿ��ٸ���Ҷ�任(FFT)����2�������������ĳ˷����㣬ʱ�临�Ӷ�O(NlogN)(��������ʱ�临�Ӷ�O(N^2))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PI 3.1415926535897932//��

#define MAX 1024//�������������������ַ��������λ��(��'\0')

int calc_N(int a_bit, int b_bit)//������ٸ���Ҷ�任��λ��(N=2^k)
{
	int N=a_bit+b_bit;
	int i=1;
	while(i<N)
	{
		i*=2;
	}
	return i;
}

class complex//����
{
public:
	complex(double x,double y);//���캯��
	double real;//ʵ��
	double image;//�鲿
	complex operator+(const complex& comp) const;
	complex operator*(const complex& comp) const;
	complex operator=(const complex& comp);//����ֵ��Ϊ�˿�������
	void operator+=(const complex& comp);//�޷���ֵ��Ϊ�˲���������
	void show();
};

complex::complex(double x=0.0,double y=0.0)
{
	real=x;
	image=y;
}

complex complex::operator+(const complex& comp) const
{
	return complex(real+comp.real, image+comp.image);
}

complex complex::operator*(const complex& comp) const
{
	return complex(real*comp.real-image*comp.image, image*comp.real+real*comp.image);
}

complex complex::operator=(const complex& comp)
{
	if(this!=&comp)
	{
		real=comp.real;
		image=comp.image;
	}
	return *this;
}

void complex::operator+=(const complex& comp)
{
	real+=comp.real;
    image+=comp.image;
}

void complex::show()//�˺���ֻ�ڵ��Դ�������в��Է���ר�ã�ֻΪ��ʾ����ǰ���ݣ����Բ�δ��0��1��-1�������������Ƚ�����ϸ������ʾ
{
	printf("%g+%gi\n",real,image);
}

complex calc_W0(int k, int N)//W0^k=e^(-j��k);��=2��/N
{
	if((4*k)%N==0)//Ϊ�˱��⸡�����������2k��+0��2k��+��/2��2k��+�У�2k��+3��/2�Ƚ��е�������(k��Z)
	{
		if((4*k/N)%4==0)
	    {
		    complex comp(1,0);
		    return comp;
	    }
	    else if((4*k/N)%4==1)
	    {
		    complex comp(0,-1);
		    return comp;
	    }
	    else if((4*k/N)%4==2)
	    {
		    complex comp(-1,0);
		    return comp;
	    }
	    else if((4*k/N)%4==3)
	    {
		    complex comp(0,1);
		    return comp;
	    }
	}
	else
	{
		complex comp(cos(-2*PI*k/N),sin(-2*PI*k/N));//ŷ����ʽ
		return comp;
	}
}

complex calc_W1(int k, int N)//W1^k=e^(j��k);��=2��/N
{
	if((4*k)%N==0)//Ϊ�˱��⸡�����������2k��+0��2k��+��/2��2k��+�У�2k��+3��/2�Ƚ��е�������(k��Z)
	{
		if((4*k/N)%4==0)
	    {
		    complex comp(1,0);
		    return comp;
	    }
	    else if((4*k/N)%4==1)
	    {
		    complex comp(0,1);
		    return comp;
	    }
	    else if((4*k/N)%4==2)
	    {
		    complex comp(-1,0);
		    return comp;
	    }
	    else if((4*k/N)%4==3)
	    {
		    complex comp(0,-1);
		    return comp;
	    }
	}
	else
	{
		complex comp(cos(2*PI*k/N),sin(2*PI*k/N));//ŷ����ʽ
	    return comp;
	}
}

double round(double r)//��������
{
    return (r>0.0)?floor(r+0.5):ceil(r-0.5);
}

int main()
{
	char a_char[MAX];//��һ������
	char b_char[MAX];//�ڶ�������
	int a_char_num=0;//��һ������λ��
	int b_char_num=0;//�ڶ�������λ��
	complex a_complex[MAX];//��һ�������Ŀ��ٸ���Ҷ�任��ֵ
	complex b_complex[MAX];//�ڶ��������Ŀ��ٸ���Ҷ�任��ֵ
	complex c_complex[2*MAX];//����Ŀ��ٸ���Ҷ�任��ֵ
	int c_int[2*MAX];//������ٸ���Ҷ��任���ֵ
	printf("�����һ��������\n");
	scanf("%s",a_char);
	printf("����ڶ���������\n");
	scanf("%s",b_char);
	a_char_num=strlen(a_char);
	b_char_num=strlen(b_char);
	char temp_char=0;
	for(int i=0;i<a_char_num/2;i++)//���ַ�������a_char����
	{
		if(a_char[i]>='0' && a_char[i]<='9')//����Ƿ�Ϊ����
		{
			temp_char=a_char[i];
		    a_char[i]=a_char[a_char_num-1-i];
		    a_char[a_char_num-1-i]=temp_char;
		}
		else
		{
			printf("������������\n");
			system("pause");
			return 0;
		}
	}
	for(int i=0;i<b_char_num/2;i++)//���ַ�������b_char����
	{
		if(b_char[i]>='0' && b_char[i]<='9')//����Ƿ�Ϊ����
		{
			temp_char=b_char[i];
		    b_char[i]=b_char[b_char_num-1-i];
		    b_char[b_char_num-1-i]=temp_char;
		}
		else
		{
			printf("������������\n");
			system("pause");
			return 0;
		}
	}
	const int N=calc_N(a_char_num,b_char_num);//������ٸ���Ҷ�任��λ��
	complex sum(0,0);
	for(int i=0;i<N;i++)//�Ե�һ���������п��ٸ���Ҷ�任
	{
		sum.real=0;
		sum.image=0;
		for(int j=0;j<a_char_num;j++)
		{
			sum+=complex(a_char[j]-'0',0)*calc_W0(j*i,N);
		}
		a_complex[i]=sum;
	}
	for(int i=0;i<N;i++)//�Եڶ����������п��ٸ���Ҷ�任
	{
		sum.real=0;
		sum.image=0;
		for(int j=0;j<b_char_num;j++)
		{
			sum+=complex(b_char[j]-'0',0)*calc_W0(j*i,N);
		}
		b_complex[i]=sum;
	}
	for(int i=0;i<N;i++)//�������Ŀ��ٸ���Ҷ�任��Ķ�Ӧֵ
	{
		c_complex[i]=a_complex[i]*b_complex[i];
	}
	for(int i=0;i<N;i++)//�Խ�����п��ٸ���Ҷ��任
	{
		sum.real=0;
		sum.image=0;
		for(int j=0;j<N;j++)
		{
			sum+=c_complex[j]*calc_W1(j*i,N);
		}
		c_int[i]=round((sum*(1.0/N)).real);
	}
	int temp_int=0;
	for(int i=0;i<N;i++)//����������������ʾ����ʽ
	{
		temp_int=c_int[i]; 
		c_int[i]=temp_int%10;
		c_int[i+1]=c_int[i+1]+temp_int/10;
	}
	printf("���Ϊ��\n");
	bool tag=false;//�Խ����ǰ��������0�Ƿ��Ѿ������ϵı�־λ
	for(int i=N-1;i>=0;i--)//�Խ���������
	{
		if(tag==true || c_int[i]!=0)
		{
			tag=true;
			printf("%d",c_int[i]);
		}
	}
	if(tag==false)//���Ϊ0
	{
		printf("0");
	}
	printf("\n");
	system("pause");
	return 0;
}