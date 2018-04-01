//此程序采用快速傅里叶变换(FFT)进行2个大整数因数的乘法运算，时间复杂度O(NlogN)(正常计算时间复杂度O(N^2))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PI 3.1415926535897932//π

#define MAX 1024//两个大整数因数读入字符串的最大位数(含'\0')

int calc_N(int a_bit, int b_bit)//计算快速傅里叶变换的位数(N=2^k)
{
	int N=a_bit+b_bit;
	int i=1;
	while(i<N)
	{
		i*=2;
	}
	return i;
}

class complex//复数
{
public:
	complex(double x,double y);//构造函数
	double real;//实部
	double image;//虚部
	complex operator+(const complex& comp) const;
	complex operator*(const complex& comp) const;
	complex operator=(const complex& comp);//返回值是为了可以连等
	void operator+=(const complex& comp);//无返回值是为了不允许连等
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

void complex::show()//此函数只在调试代码过程中测试方便专用，只为显示出当前数据，所以并未对0，1，-1，正数，负数等进行详细分类显示
{
	printf("%g+%gi\n",real,image);
}

complex calc_W0(int k, int N)//W0^k=e^(-jωk);ω=2π/N
{
	if((4*k)%N==0)//为了避免浮点型运算的误差，2kπ+0，2kπ+π/2，2kπ+π，2kπ+3π/2等进行单独处理(k∈Z)
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
		complex comp(cos(-2*PI*k/N),sin(-2*PI*k/N));//欧拉公式
		return comp;
	}
}

complex calc_W1(int k, int N)//W1^k=e^(jωk);ω=2π/N
{
	if((4*k)%N==0)//为了避免浮点型运算的误差，2kπ+0，2kπ+π/2，2kπ+π，2kπ+3π/2等进行单独处理(k∈Z)
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
		complex comp(cos(2*PI*k/N),sin(2*PI*k/N));//欧拉公式
	    return comp;
	}
}

double round(double r)//四舍五入
{
    return (r>0.0)?floor(r+0.5):ceil(r-0.5);
}

int main()
{
	char a_char[MAX];//第一个因数
	char b_char[MAX];//第二个因数
	int a_char_num=0;//第一个因数位数
	int b_char_num=0;//第二个因数位数
	complex a_complex[MAX];//第一个因数的快速傅里叶变换的值
	complex b_complex[MAX];//第二个因数的快速傅里叶变换的值
	complex c_complex[2*MAX];//结果的快速傅里叶变换的值
	int c_int[2*MAX];//结果快速傅里叶逆变换后的值
	printf("输入第一个因数：\n");
	scanf("%s",a_char);
	printf("输入第二个因数：\n");
	scanf("%s",b_char);
	a_char_num=strlen(a_char);
	b_char_num=strlen(b_char);
	char temp_char=0;
	for(int i=0;i<a_char_num/2;i++)//对字符串数组a_char逆序
	{
		if(a_char[i]>='0' && a_char[i]<='9')//检测是否为数字
		{
			temp_char=a_char[i];
		    a_char[i]=a_char[a_char_num-1-i];
		    a_char[a_char_num-1-i]=temp_char;
		}
		else
		{
			printf("输入数据有误！\n");
			system("pause");
			return 0;
		}
	}
	for(int i=0;i<b_char_num/2;i++)//对字符串数组b_char逆序
	{
		if(b_char[i]>='0' && b_char[i]<='9')//检测是否为数字
		{
			temp_char=b_char[i];
		    b_char[i]=b_char[b_char_num-1-i];
		    b_char[b_char_num-1-i]=temp_char;
		}
		else
		{
			printf("输入数据有误！\n");
			system("pause");
			return 0;
		}
	}
	const int N=calc_N(a_char_num,b_char_num);//计算快速傅里叶变换的位数
	complex sum(0,0);
	for(int i=0;i<N;i++)//对第一个因数进行快速傅里叶变换
	{
		sum.real=0;
		sum.image=0;
		for(int j=0;j<a_char_num;j++)
		{
			sum+=complex(a_char[j]-'0',0)*calc_W0(j*i,N);
		}
		a_complex[i]=sum;
	}
	for(int i=0;i<N;i++)//对第二个因数进行快速傅里叶变换
	{
		sum.real=0;
		sum.image=0;
		for(int j=0;j<b_char_num;j++)
		{
			sum+=complex(b_char[j]-'0',0)*calc_W0(j*i,N);
		}
		b_complex[i]=sum;
	}
	for(int i=0;i<N;i++)//计算结果的快速傅里叶变换后的对应值
	{
		c_complex[i]=a_complex[i]*b_complex[i];
	}
	for(int i=0;i<N;i++)//对结果进行快速傅里叶逆变换
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
	for(int i=0;i<N;i++)//将结果整理成正常表示的形式
	{
		temp_int=c_int[i]; 
		c_int[i]=temp_int%10;
		c_int[i+1]=c_int[i+1]+temp_int/10;
	}
	printf("结果为：\n");
	bool tag=false;//对结果最前方连续个0是否已经输出完毕的标志位
	for(int i=N-1;i>=0;i--)//对结果进行输出
	{
		if(tag==true || c_int[i]!=0)
		{
			tag=true;
			printf("%d",c_int[i]);
		}
	}
	if(tag==false)//结果为0
	{
		printf("0");
	}
	printf("\n");
	system("pause");
	return 0;
}