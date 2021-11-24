// CSCI570FinalProject.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>    // std::min
#include <vector>
using namespace std;

//Global Variable
const int GapPenalty = 30;
int const MissPenalityArray[4][4] = {
    {0, 110, 48, 94},
    {110, 0, 118, 48},
    {48, 118, 0, 110},
    {94, 48, 110, 0}
};
bool is_number(const string& s)
{
    string::const_iterator it = s.begin();
    while (it != s.end() && isdigit(*it)) it++;
    return !s.empty() && it == s.end();
}

string* generateInput()
{
    string gene1, gene2;
    string input;
    string* result = new string[2];
    int index;
    ifstream nameFileout;
    nameFileout.open("input.txt");
    while (getline(nameFileout, input))
    {
        if (is_number(input))
        {
            index = stoi(input) + 1;
            if (gene2.empty())
            {
                gene1.insert(index, gene1);
            }
            else {
                gene2.insert(index, gene2);
            }
        }
        else {
            if (gene1.empty())
            {
                gene1 = input;
            }
            else 
            {
                gene2 = input;
            }
        }
    }
    result[0] = gene1;
    result[1] = gene2;
    return result;
}

int calculateMissPenality(char x, char y)
{
//    return 3;
    int index1, index2;
    switch (x) {
    case 'A':
        index1 = 0;
        break;
    case 'C':
        index1 = 1;
        break;
    case 'G':
        index1 = 2;
        break;
    case 'T':
        index1 = 3;
        break;
    default: index1 = 0;
        break;
    };
    switch (y) {
    case 'A':
        index2 = 0;
        break;
    case 'C':
        index2 = 1;
        break;
    case 'G':
        index2 = 2;
        break;
    case 'T':
        index2 = 3;
        break;
    default: index2 = 0;
        break;
    };
    return MissPenalityArray[index1][index2];
}

int* DP_SeqAlignment(string gene1, string gene2)
{
    const size_t m = gene1.length();
    const size_t n = gene2.length();
    int* dp = new int[(m + 1) * (n + 1)];

    dp[0] = 0;
    for (int i = 0; i <= m; i++)                            //initialize first row to max gap penality
    {
        dp[i] = i * GapPenalty;
    }

    for (int j = 0; j <= n; j++)                            //initialize first column to max gap penality
    {
        dp[j * (m + 1)] = j * GapPenalty;
    }

    for (int i = 1; i <= m; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            if (gene1[i - 1] == gene2[j - 1])
            {
                dp[(m + 1) * j + i] = dp[(m + 1) * (j - 1) + (i - 1)];
            }
            else
            {
                dp[(m + 1) * j + i] = min(dp[(m + 1) * (j - 1) + (i - 1)] + calculateMissPenality(gene1[i - 1], gene2[j - 1]),
                                          min(dp[(m + 1) * j + (i - 1)] + GapPenalty,
                                          dp[(m + 1) * (j - 1) + i] + GapPenalty));
            }
        }
    }
    return dp;
}

void constructAlignment(string gene1, string gene2, int* dp)
{
    const size_t m = gene1.length();
    const size_t n = gene2.length();
    int l = n + m;
    int i = m;
    int j = n;
    int xpos = l;
    int ypos = l;
    int* xAns = new int[l + 1];
    int* yAns = new int[l + 1];

    while (!(i == 0 || j == 0))
    {
        if (gene1[i - 1] == gene2[j - 1])
        {
            xAns[xpos--] = (int)gene1[i - 1];
            yAns[ypos--] = (int)gene2[j - 1];
            i--; j--;
        }
        else if ((dp[(m + 1) * (j - 1) + (i - 1)] + calculateMissPenality(gene1[i - 1], gene2[j - 1])) == dp[(m + 1) * j + i])
        {
            xAns[xpos--] = (int)gene1[i - 1];
            yAns[ypos--] = (int)gene2[j - 1];
            i--; j--;
        }
        else if (dp[(m + 1) * j + (i - 1)] + GapPenalty == dp[(m + 1) * j + i])
        {
            xAns[xpos--] = (int)gene1[i - 1];
            yAns[ypos--] = (int)'_';
            i--;
        }
        else if ((dp[(m + 1) * (j - 1) + i] + GapPenalty) == dp[(m + 1) * j + i])
        {
            xAns[xpos--] = (int)'_';
            yAns[ypos--] = (int)gene2[j - 1];
            j--;
        }
    }
    while (xpos > 0)
    {
        if (i > 0)
            xAns[xpos--] = (int)gene1[--i];
        else
            xAns[xpos--] = (int)'_';
    }
    while (ypos > 0)
    {
        if (j > 0)
            yAns[ypos--] = (int)gene2[--j];
        else
            yAns[ypos--] = (int)'_';
    }

    int id = 1;
    for (i = l; i >= 1; i--)
    {
        if ((char)yAns[i] == '_' && (char)xAns[i] == '_')
        {
            id = i + 1;
            break;
        }
    }
    for (int i = 0; i <= m; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            cout << dp[j * (m + 1) + i] << ' ';
        }
        cout << endl;
    }
    cout << "Minimum Penalty in aligning the genes = ";
    cout << dp[((m + 1) * (n)) + m] << "\n";

    for (i = id; i <= l; i++)
    {
        cout << (char)xAns[i];
    }
    cout << "\n";
    for (i = id; i <= l; i++)
    {
        cout << (char)yAns[i];
    }
}

int* MmeoryEfficientAlignment(string gene1, string gene2)
{
    const size_t m = gene1.length();
    const size_t n = gene2.length();
    int* dp;
    if (m <= 2 || n <= 2)
    {
        dp = DP_SeqAlignment(gene1, gene2);
        return dp;
    }
    int MidIndex = gene1.length() / 2;
    string gene1_left = gene1.substr(0, MidIndex);
    int* dp_left = DP_SeqAlignment(gene1_left, gene2);

    string gene1_right = gene1.substr(MidIndex, m-MidIndex);
    string gene1_right_reverse = string(gene1_right.rbegin(), gene1_right.rend());
    string gene2_reverse = string(gene2.rbegin(), gene2.rend());
    int* dp_right = DP_SeqAlignment(gene1_right_reverse, gene2_reverse);
    
    
}
int main()
{
    string* input = generateInput();
    cout << input[0] << endl;
    cout << input[1] << endl;
    int* dp = DP_SeqAlignment(input[0], input[1]);
    constructAlignment(input[0], input[1], dp);
}