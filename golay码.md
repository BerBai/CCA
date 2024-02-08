
### 码字结构

**Golay (23,12)**

表示该编码共有 23 个位、12 个信息位和 23-12=11 个校验位，由于每个码字的长度为 23 位，因此可能的二进制数值有 $2^{23}$ 个或 8,388,608 个。但是，由于每个 12 位信息字段只有一组相应的 11 个校验位，因此只有 $2^{12}$ 个或 4096 个有效的 Golay 码字。

![](Pasted%20image%2020240204152458.png)

golay 码字并不是像沿着数字线数的那样均匀分布，比如每隔几个刻度。而是根据它们之间的汉明距离进行间隔。事实证明，每个戈莱码字都有七个或七个以上的比特相互不同。Golay (23,12)编码的最小距离 $d$ 为 7，在任何模式下最多可检测并纠正 $(d-1)/2=3$ 位 bit 错误。

**Golay (24,12)**

为了增强戈莱码的能力，通常会添加一个奇偶校验位，从而产生一个简洁的三字节码元，称为扩展 golay 码。

![](Pasted%20image%2020240204154923.png)

![](Pasted%20image%2020240204155742.png)

上图说明了码字校正过程中发生的情况。我们看到两个有效的 golay码字，即 A 和 B。Y 轴代表远离有效码字时产生的比特错误数。将校正算法想象成获取损坏的输入码字（红色标记），然后沿着斜坡滑向最近的正确码字。如果该编解码是编解码 B 的损坏版本，我们就走运了。如果它是 A 的严重破坏版本，校正算法就会欺骗我们，返回 B 作为校正后的编解码。我们可以看到，编码的最小距离 d 控制着可以容忍的损坏程度。当然，这只是一个简化的解释。

二元golay 码字具有一些非常有趣的特性：
- 循环不变性。如果将一个 23 位的 golay 码字循环移位任意位数，其结果也是一个有效的 golay码字。
- 反转性。如果将 23 位 golay 码字反转，其结果也是有效的 golay码字。
- 最小汉明距离。任何两个 golay (23,12)码元之间的距离总是 7 位或更多位。任何两个 golay (24,12) 码字之间的距离总是八个或更多比特。
- 纠错。纠错算法可检测并纠正每个码字最多三个比特的错误。

三元 Golay 码，它的操作对象是三元而非二元数字。三元 Golay 码将每6个三元符号分为一组，编码生成5个冗余校验三元符号。这样由11个三元符号组成的三元 Golay 码码字可以纠正2个错误。


### 代码实现

#### 编码

Golay 代码的编码方式使用模 2 除法。但是，每个码字只有 12 个信息位，因此您必须将数据分解为 12 位块，并将每个块编码为一个码字。Golay 代码的特征多项式是：

- $X^{11}+ X^9+ X^7+ X^6+ X^5+ X + 1$，系数 `0xAE3`。
- $X^{11}+ X^{10}+ X^6+ X^5+ X^4+ X^2+ 1$，系数 `0xC75`。

有两个多项式，只需要使用其中之一，因为它们都具有与上述相同的属性，尽管它们生成不同的校验位。我们将使用第一个多项式 AE3h 的 X 项的系数。

请注意，`0xAE3` 是 `0xC75` 的位反转版本，是代码循环性质的副产品。

编码的

```c
/******************* main.cpp *******************/  
  
#include <iostream>  
#include <vector>  
#include <math.h>  
#include <time.h>  
#include <string>  
#include "golay.h"  
#include <ctime>  
  
using namespace std;  
  
int main() {  
    system("chcp 65001");  
  
    matrix gen,gen_t,par,par_t,encoded_message,syndrome,  
            sB,s_add_colB,sB_add_colB_T,original_message;  
    matrix message(1,12);  
    matrix temp1(1,12);  
    matrix noise(1,24);  
    matrix error(1,24);  
    matrix s_add_colB_full(12,13);  
    matrix sB_add_colB_T_full(12,13);  
    matrix temp(12,1);  
  
    int i,j,k;  
    vector<int> zeros(24);  
  
    clock_t start = clock();  
  
    gen = golay_generator();   // 生成 Golay 24 码的生成矩阵  
    cout << "生成 Golay 24 码的生成矩阵（G）：" << endl;  
    print(gen);  
  
    gen_t = transpose(gen); // 生成矩阵 G 的转置  
    cout << endl << "生成矩阵转置（G^T）：" << endl;  
    print(gen_t);  
  
    par = generate_B();        // Golay 24 码的奇偶校验矩阵  
    cout << endl << "Golay 24 码的奇偶校验矩阵（B）：" << endl;  
    print(par);  
  
    par_t = transpose(par); // 奇偶校验矩阵 B 的转置  
    cout << endl << "奇偶校验矩阵转置（B^T）：" << endl;  
    print(par_t);  
  
    message = generate_message();   // 生成随机消息  
    cout << endl << "消息（m）：   ";  
    print(message);  
  
    encoded_message = encode(message,gen);  // 编码消息  
    cout << endl << "编码消息（mG）：  ";  
    print(encoded_message);  
  
    noise = noisy_channel(encoded_message); // 将消息发送到“有噪声的信道”  
    cout << endl << "接收到的消息（r = c + e）：  ";  
    print(noise);  
  
    syndrome = multiply(noise, gen_t);  // 计算 s = rG^T    cout << endl << "综合（s = rG^T）：   ";  
    print(syndrome);  
    cout << "综合的权重是 " << weight(syndrome) << endl;  
  
    sB = multiply(syndrome,par);    // 计算 sB    cout << endl << "sB：    ";  
    print(sB);  
    cout  << "sB 的权重是 " << weight(sB);  
  
    /* 计算行向量 s + (c_j)^T 并计算行向量的权重 */    for(i = 0 ; i < 12 ; i++) {  
        for(j = 0 ; j < 12 ; j++) {  
            temp.m[j][0] = par.m[j][i];  
        }  
  
        temp1 = transpose(temp);  
        s_add_colB = addition(syndrome,temp1);  
  
        for(k = 0 ; k < 13 ; k++) {  
            if(k < 12) {  
                s_add_colB_full.m[i][k] = s_add_colB.m[0][k];  
            } else {  
                s_add_colB_full.m[i][k] = weight(s_add_colB);  
            }  
        }  
    }  
  
    cout << "\n\n" << "完整的 s + (c_j)^T（最后一列的值是行的权重）：" << endl;  
    print(s_add_colB_full);  
  
    /* 计算行向量 sB + (b_j)^T 并计算行向量的权重 */    for(i = 0 ; i < 12 ; i++) {  
        for(j = 0 ; j < 12 ; j++) {  
            temp.m[j][0] = par_t.m[j][i];  
        }  
  
        temp1 = transpose(temp);  
        sB_add_colB_T = addition(sB,temp1);  
  
        for(k = 0 ; k < 13 ; k++) {  
            if(k < 12) {  
                sB_add_colB_T_full.m[i][k] = sB_add_colB_T.m[0][k];  
            } else {  
                sB_add_colB_T_full.m[i][k] = weight(sB_add_colB_T);  
            }  
        }  
    }  
  
    cout << "\n\n" << "完整的 sB + (b_j)^T（最后一列的值是行的权重）：" << endl;  
    print(sB_add_colB_T_full);  
  
    /***** 计算错误向量 e *****/    if(weight(syndrome) <= 3) {  
        for(i = 0; i < syndrome.m[0].size() ; i++) {  
            if(syndrome.m[0][i] == 1) {  
                error.m[0][i] = 1;  
            }  
        }  
    } else if(weight(sB) <= 3) {  
        for(i = 0; i < sB.m[0].size() ; i++) {  
            if(sB.m[0][i] == 1) {  
                error.m[0][i+12] = 1;  
            }  
        }  
    }  
  
    if(error.m[0] == zeros) {  
        for(i = 0; i < 12 ; i++) {  
            if(s_add_colB_full.m[i][12] <= 2) {  
                error.m[0][i+12] = 1;  
  
                for(j = 0 ; j< 12 ; j++) {  
                    if(s_add_colB_full.m[i][j] ==1) {  
                        error.m[0][j] = 1;  
                    }  
                }  
            }  
        }  
    }  
  
    if(error.m[0] == zeros) {  
        for(i = 0; i < 12 ; i++) {  
            if(sB_add_colB_T_full.m[i][12] <= 2) {  
                error.m[0][i] = 1;  
  
                for(j = 0 ; j< 12 ; j++) {  
                    if(sB_add_colB_T_full.m[i][j] ==1) {  
                        error.m[0][j+12] = 1;  
                    }  
                }  
            }  
        }  
    }  
  
    cout << endl <<"错误向量是：    ";  
    print(error);  
  
    original_message =addition(noise,error);    // 计算原始消息 c = r + e  
    /***** 检查解码是否正确 *****/    if(original_message.m == encoded_message.m) {  
        cout << endl << "传输和纠错成功。" << endl;  
    } else {  
        cout << "传输或纠错出现错误。" << endl;  
    }  
  
    cout << endl << "计算出的原始发送消息：   ";  
    print(original_message);  
    cout << "         原始发送消息：   ";  
    print(encoded_message);  
  
  
    cout << endl << "计算出的原始消息 m ：   ";  
    for(i = 0 ; i < 12 ; i++) {  
        cout << original_message.m[0][i] << "  ";  
    }  
    cout << endl << "         原始消息 m ：   ";  
    for(i = 0 ; i < 12 ; i++) {  
        cout << message.m[0][i] << "  ";  
    }  
  
    clock_t finish = clock();  
  
    cout << "\n\n" << "运行时间： " << (double)(finish-start)*1000/CLOCKS_PER_SEC<< " ms" << endl;  
  
    return 0;  
}

```

```cpp
/*********************golay.h*********************/  
  
#include <iostream>  
#include <vector>  
#include <math.h>  
#include <time.h>  
#include <string>  
  
using namespace std;  
  
class matrix  
{  
public:  
    int row;  
    int column;  
    vector < vector <int> > m;  
  
    matrix() {};  
    matrix(int r, int c)  
    {  
        row = r;  
        column = c;  
  
        for (int i = 0; i < r; i++)  
        {  
            m.push_back(vector<int>());  
        }  
  
        for (int j = 0; j < c; j++)  
        {  
            for (int i = 0; i < m.size(); i++)  
            {  
                m[i].push_back(0);  
            }  
        }  
    }  
};  
  
// 矩阵乘法  
matrix multiply(matrix a, matrix b);  
// 打印矩阵  
void print(matrix a);  
// 转置矩阵  
matrix transpose(matrix a);  
// 生成Golay码的生成矩阵  
matrix golay_generator();  
// 生成B矩阵  
matrix generate_B();  
// 生成单位矩阵  
matrix generate_I(int a);  
// 矩阵加法  
matrix addition(matrix a, matrix b);  
// 生成消息  
matrix generate_message();  
// 编码  
matrix encode(matrix m);  
// 产生噪声  
matrix noisy_channel(matrix m);  
// 权重  
int weight(matrix a);  
  
// 矩阵乘法  
matrix multiply(matrix a, matrix b)  
{  
    int i,j,k,l;  
    int tem = 0;  
    int index = 0;  
    matrix temp(a.row , b.column);  
  
    if (a.column != b.row)  
    {  
        cout << "错误：无法相乘这两个矩阵";  
    }  
    else  
    {  
        for(i = 0 ; i < a.row ; i++)  
        {  
            for(j = 0; j < b.column ; j++)  
            {  
                for(k = 0 ; k < a.column ; k++)  
                {  
                    tem = a.m[i][k] * b.m[k][j];  
                    index = index + tem;  
                    tem = 0;  
                }  
  
                temp.m[i][j] = index % 2;  
                index = 0;  
            }  
        }  
    }  
  
    return temp;  
}  
  
// 矩阵加法  
matrix addition(matrix a, matrix b)  
{  
    matrix temp(a.row, a.column);  
    int i,j;  
  
    if(a.row != b.row || a.column != b.column)  
    {  
        cout << "错误：无法相加这两个矩阵";  
    }  
    else  
    {  
        for(i = 0 ; i < a.row ; i++)  
        {  
            for(j = 0 ; j < a.column ; j++)  
            {  
                temp.m[i][j] = (a.m[i][j] + b.m[i][j]) % 2;  
            }  
        }  
    }  
  
    return temp;  
}  
  
// 矩阵转置  
matrix transpose(matrix a)  
{  
    matrix temp(a.column , a.row);  
    int i, j;  
  
    for(i = 0; i < a.row ; i++)  
    {  
        for(j = 0 ; j < a.column ; j++)  
        {  
            temp.m[j][i] = a.m[i][j];  
        }  
    }  
  
    return temp;  
}  
  
// 打印矩阵  
void print(matrix a)  
{  
    int i,j;  
  
    for (i = 0 ; i < a.row ; i++)  
    {  
        for (j = 0 ; j < a.column ; j++)  
        {  
  
            cout << a.m[i][j] << "  ";  
        }  
  
        cout << endl;  
    }  
}  
  
// 生成Golay码的生成矩阵  
matrix golay_generator()  
{  
    matrix generator(12,24);  
    int i, j, k;  
    int check = 0;  
    int count = 1;  
    vector<int> squares(1);  
  
    for(i = 0 ; i < 12 ; i++)  
    {  
        generator.m[i][i] = 1;  
    }  
  
    for(i = 0 ; i < 11 ; i++)  
    {  
        generator.m[i][12] = 1;  
    }  
  
    for(i = 13 ; i < 24 ; i++)  
    {  
        generator.m[11][i] = 1;  
    }  
  
    for(i = 0 ; i < 11 ; i++)  
    {  
        j = (i*i) % 11;  
  
        for(k = 0; k < squares.size() ; k++)  
        {  
            if(j == squares[k])  
            {  
                check = 1;  
            }  
  
        }  
  
        if(check == 0)  
        {  
            squares.push_back(j);  
        }  
  
        check = 0;  
    }  
  
    for(i = 0 ; i < squares.size() ; i++)  
    {  
        squares[i] = squares[i] + 13;  
    }  
  
    for(i = 13 ; i < 24 ; i++)  
    {  
        for(j = 0 ; j < squares.size() ; j++)  
        {  
            if(squares[j] == i)  
            {  
                generator.m[0][i] = 1;  
            }  
        }  
    }  
  
    while(count != 11)  
    {  
        for(i = 0 ; i < squares.size() ; i++)  
        {  
            squares[i] = ((squares[i] - 12) % 11) + 13;  
        }  
  
        for(i = 13 ; i < 24 ; i++)  
        {  
            for(j = 0 ; j < squares.size() ; j++)  
            {  
                if(squares[j] == i)  
                {  
                    generator.m[count][i] = 1;  
                }  
            }  
        }  
  
        count++;  
    }  
  
    return generator;  
}  
  
// 生成B矩阵  
matrix generate_B()  
{  
    matrix generator(12,12);  
    int i, j, k;  
    int count = 1;  
    int check = 0;  
    vector<int> squares;  
  
    for(i = 0 ; i < 11 ; i++)  
    {  
        generator.m[i][0] = 1;  
    }  
  
    for(i = 1 ; i < 12 ; i++)  
    {  
        generator.m[11][i] = 1;  
    }  
  
    for(i = 0 ; i < 11 ; i++)  
    {  
        j = (i*i) % 11;  
  
        for(k = 0; k < squares.size() ; k++)  
        {  
            if(j == squares[k])  
            {  
                check = 1;  
            }  
  
        }  
  
        if(check == 0)  
        {  
            squares.push_back(j);  
        }  
  
        check = 0;  
    }  
  
    for(i = 0 ; i < squares.size() ; i++)  
    {  
        squares[i] = squares[i] + 1;  
    }  
  
    for(i = 0 ; i < 12 ; i++)  
    {  
        for(j = 0 ; j < squares.size() ; j++)  
        {  
            if(squares[j] == i)  
            {  
                generator.m[0][i] = 1;  
            }  
        }  
    }  
  
    while(count != 11)  
    {  
        for(i = 0 ; i < squares.size() ; i++)  
        {  
            squares[i] = (squares[i] % 11) + 1;  
  
  
        }  
  
        for(i = 1 ; i < 12 ; i++)  
        {  
            for(j = 0 ; j < squares.size() ; j++)  
            {  
                if(squares[j] == i)  
                {  
                    generator.m[count][i] = 1;  
                }  
            }  
        }  
  
        count++;  
    }  
  
    return generator;  
}  
  
// 生成单位矩阵  
matrix generate_I(int a)  
{  
    matrix temp(a,a);  
    int i;  
  
    for(i = 0 ; i < a ; i++)  
    {  
        temp.m[i][i] = 1;  
    }  
  
    return temp;  
}  
  
// 生成消息  
matrix generate_message()  
{  
    matrix message(1,12);  
    int i,temp;  
  
    srand (time(NULL));  
  
    for(i = 0 ; i < 12 ; i++)  
    {  
        temp = rand()%2;  
        message.m[0][i] = temp;  
    }  
  
    return message;  
}  
  
// 编码  
matrix encode(matrix m, matrix G)  
{  
    matrix temp = multiply(m,G);  
  
    return temp;  
}  
  
// 产生噪声  
matrix noisy_channel(matrix m)  
{  
    int i = 0;  
    int r;  
    int errors = 0;  
    matrix temp = m;  
  
    srand (time(NULL));  
  
    while(errors != 3 && i != temp.m[0].size())  
    {  
  
        r = rand()%100;  
  
        if(r >= 25 && r <= 35)  
        {  
            temp.m[0][i] = (temp.m[0][i] + 1) % 2;  
  
            errors++;  
        }  
  
        i++;  
    }  
  
    cout << endl << "传输过程中的错误数量为：" << errors << endl;  
  
    return temp;  
}  
  
// 权重  
int weight(matrix a)  
{  
    int temp = 0;  
    int i;  
  
    for(i = 0 ; i < a.m[0].size() ; i++)  
    {  
        if(a.m[0][i] == 1)  
        {  
            temp++;  
        }  
    }  
  
    return temp;  
}
```