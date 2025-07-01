#include <iostream>
#include <vector>
#include <limits>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <map>
#include <bits/stdc++.h>

#define RESET "\033[0m"
#define RED "\033[31m"
#define BLUE "\033[34m"
#define GREEN "\033[32m"
#define YELLOW "\033[33m"
#define VIOLET "\033[35m"

using namespace std;

const int M = 100100;
bool UnboundedSol = false;
bool Nosol = false;
bool alternateSol = false;
bool degeneracy = false;
bool degeneracyType1 = false;
bool infeasiblesol = false;
vector<map<int, double>> saveAlternatesol;
bool fhase1Complte = false;
bool minTypeOptimalSol = false;

const double INF = numeric_limits<double>::infinity();
const double EPS = 1e-9;
vector<int> artificialval;

int numberOfVariabe;
int numberOfEquation;

double findzj(vector<double> c, vector<vector<double>> a, vector<int> bi, int ith)
{
    double result = 0;
    for (int i = 0; i < bi.size(); i++)
    {
        result += c[bi[i]] * a[i][ith];
    }
    return result;
}

void operation(int pivoti, int pivotj, vector<vector<double>> &a, vector<double> &b)
{
    // cout<<"\ntest 1   ";
    double pivot = a[pivoti][pivotj];
    for (int i = 0; i < a[pivoti].size(); i++)
    {
        a[pivoti][i] /= pivot;
        // cout<<a[pivoti][i]<<"  ";
    }

    b[pivoti] /= pivot;

    for (int i = 0; i < a.size(); i++)
    {
        if (i != pivoti)
        {
            double tem = a[i][pivotj];
            for (int j = 0; j < a[i].size(); j++)
            {
                a[i][j] = a[i][j] - (a[pivoti][j] * tem);
            }

            b[i] = b[i] - b[pivoti] * tem;
        }
    }
}

vector<int> checkDegeracy(double minxbyr, vector<double> xbyr)
{
    vector<int> repeat;
    for (int i = 0; i < xbyr.size(); i++)
    {
        if (xbyr[i] == minxbyr)
        {
            repeat.push_back(i);
        }
        // cout << "\n"
        //      << xbyr[i] << endl;
    }

    return repeat;
}

void printTopLine(vector<double> c)
{
    cout << "    |    |    " << " Cj\t";
    for (int i = 0; i < c.size(); i++)
    {
        if (abs(c[i]) >= 100000)
        {
            cout << " -M\t";
        }
        else if (c[i] < 0 || c[i] / 10 > 0)
            cout << c[i] << '\t';
        else
            cout << c[i] << '\t';
    }
    cout << "\n";
}

void printLowerLine(vector<double> zjcj)
{
    cout << "    |     Zj-Cj\t";
    // for (int i = 0; i < zjcj.size(); i++)
    // {
    //     cout << zjcj[i] << ' ';
    // }
    for (int i = 0; i < zjcj.size(); i++)
    {
        // Alternating color logic
        if (i % 2 == 0)
        {
            // Print in blue for even indices
            cout << BLUE << zjcj[i] << '\t' << RESET;
        }
        else
        {
            // Print in green for odd indices
            cout << GREEN << zjcj[i] << '\t' << RESET;
        }
    }
    cout << '\n';
}

void printMidLeft(const vector<int> &basicindex, const vector<double> &c, const vector<vector<double>> &a, const vector<double> &b, const vector<int> &repeat, int checkxi, int pivotj, vector<double> &xbyrvec, double &minxbyr, int &pivoti, map<int, double> &basicsolution)
{
    for (int i = 0; i < basicindex.size(); i++)
    {
        cout << ' ' << c[basicindex[i]] << "  ";
        cout << " a" << basicindex[i] + 1 << " ";
        cout << " x" << basicindex[i] + 1 << " ";
        basicsolution[basicindex[i]] = b[i];
        cout << " " << b[i] << " ";
        // for (int j = 0; j < c.size(); j++)
        // {
        //     if (a[i][j] < 0)
        //         cout << " " << a[i][j] << " ";
        //     else
        //         cout << "  " << a[i][j] << " ";
        // }

        for (int j = 0; j < a[i].size(); j++)
        {
            // Alternating color logic
            if ((j) % 2 == 0)
            {
                // Print in blue for even positions
                if (a[i][j] < 0)
                    cout << BLUE << " " << a[i][j] << " " << RESET;
                else
                    cout << BLUE << "  " << a[i][j] << " " << RESET;
            }
            else
            {
                // Print in green for odd positions
                if (a[i][j] < 0)
                    cout << GREEN << " " << a[i][j] << " " << RESET;
                else
                    cout << GREEN << "  " << a[i][j] << " " << RESET;
            }
        }
        // cout << endl;

        if (find(repeat.begin(), repeat.end(), i) != repeat.end())
        {
            double ratio = a[i][checkxi] / a[i][pivotj];
            cout << YELLOW << " " << ratio << RESET;
            xbyrvec.push_back(ratio);

            if (a[i][pivotj] && ratio >= 0 && minxbyr > ratio)
            {
                minxbyr = ratio;
                pivoti = i;
            }
        }
        cout << '\n';
    }
}

void printMidRight(vector<double> c)
{
    cout << " CB | B  | XB " << " b\t";
    for (int i = 0; i < c.size(); i++)
    {
        cout << "X" << i + 1 << '\t';
    }
    cout << "Xb/yr\n";
}

vector<int> checkAlternateSol(vector<double> zjcj, vector<int> basicindex)
{
    vector<int> zero;
    for (int i = 0; i < zjcj.size(); i++)
    {
        if (zjcj[i] == 0)
            zero.push_back(i);
    }
    return zero;
}

int checkNonbasicShowzero(vector<int> zeroIndex, vector<int> basicindex, vector<vector<double>> a, int pivoti)
{
    cout << "zero : ";
    for (auto i : zeroIndex)
    {
        cout << i << ' ';
    }
    cout << "\nbasic : ";
    for (auto i : basicindex)
    {
        cout << i << ' ';
    }
    cout << "\n";
    for (auto i : zeroIndex)
    {
        if (find(basicindex.begin(), basicindex.end(), i) == basicindex.end())
        {
            cout << "Non basic solution that show {zj-cj == 0} are : " << i + 1 << endl;
            if (a[pivoti][i] != 0)
                return i;
        }
    }
    return -1;
}

void printtable(vector<vector<double>> &a, vector<double> &b, vector<double> &c, vector<int> &basicindex)
{
    // cout << fixed << setprecision(4);

    bool isnosol = true;

    for (int i = 0; i < 20; i++)
    {
        map<int, double> basicsolution;

        printTopLine(c);

        printMidRight(c);

        // double zj = findzj(c,a,basicindex,0);

        vector<double> zjcj;
        int pivoti = 0, pivotj = 0;
        double minzjcj = 100100;
        double minxbyr = 100100;

        for (int i = 0; i < c.size(); i++)
        {
            double val = findzj(c, a, basicindex, i) - c[i];
            zjcj.push_back(val);
            if (minzjcj > val)
            {
                pivotj = i;
                minzjcj = val;
            }
        }

        // checkDegeneracy()

        vector<double> xbyrvec;

        basicsolution.clear();

        for (int i = 0; i < basicindex.size(); i++)
        {

            if (abs(c[basicindex[i]]) == M)
            {
                cout << " -M ";
            }
            else
            {
                cout << ' ' << c[basicindex[i]] << "   ";
            }
            cout << " a" << basicindex[i] + 1 << " ";
            cout << "  x" << basicindex[i] + 1 << " ";
            cout << " " << b[i] << "\t";
            basicsolution[basicindex[i]] = b[i];

            // cout<<" . "<<"("<<basicindex[i]<<","<<basicsolution[basicindex[i]]<<") "<<b[i];

            // for (int j = 0; j < c.size(); j++)
            // {
            //     if (a[i][j] < 0)
            //         cout << " " << a[i][j] << " ";
            //     else
            //         cout << "  " << a[i][j] << " ";
            // }
            for (int j = 0; j < a[i].size(); j++)
            {
                // Alternating color logic
                if ((j) % 2 == 0)
                {
                    // Print in blue for even positions
                    if (a[i][j] < 0)
                        cout << BLUE << a[i][j] << "\t" << RESET;
                    else
                        cout << BLUE << a[i][j] << "\t" << RESET;
                }
                else
                {
                    // Print in green for odd positions
                    if (a[i][j] < 0)
                        cout << GREEN << a[i][j] << "\t" << RESET;
                    else
                        cout << GREEN << a[i][j] << "\t" << RESET;
                }
            }
            // cout << endl;

            cout << YELLOW << b[i] / a[i][pivotj] << "\t" << RESET;
            xbyrvec.push_back(b[i] / a[i][pivotj]);
            if (a[i][pivotj] > 0 && (b[i] / a[i][pivotj]) >= 0 && minxbyr > (b[i] / a[i][pivotj]))
            {
                minxbyr = b[i] / a[i][pivotj];
                pivoti = i;
            }
            cout << '\n';
        }

        // for(auto i : basicsolution){
        //     cout<<"("<<i.first<<","<<i.second<<") ";
        // }cout<<endl;

        if (minxbyr == 100100)
        {
            cout << RED << "\n\nUnbounded solution\n\n"
                 << RESET;
            UnboundedSol = true;
            return;
        }

        cout << "    |     Zj-Cj\t\t";
        // for (int i = 0; i < zjcj.size(); i++)
        // {
        //     cout << zjcj[i] << ' ';
        // }

        for (int i = 0; i < zjcj.size(); i++)
        {
            // Alternating color logic
            if (i % 2 == 0)
            {
                // Print in blue for even indices
                cout << BLUE << zjcj[i] << '\t' << RESET;
            }
            else
            {
                // Print in green for odd indices
                cout << GREEN << zjcj[i] << '\t' << RESET;
            }
        }
        // vector<int
        vector<int> repeat = checkDegeracy(minxbyr, xbyrvec);

        if (repeat.size() > 1)
        {
            int checkxi = numberOfVariabe;

            vector<int> oldrepeat = repeat;

            for (; (repeat.size() > 1); checkxi++)
            {

                repeat = oldrepeat;
                cout << '\n'
                     << "\nSolution show Degeneracy\n\n";
                degeneracy = true;

                printTopLine(c);

                cout << " CB | B  | XB " << " b  ";
                for (int i = 0; i < c.size(); i++)
                {
                    cout << " X" << i + 1 << ' ';
                }
                cout << " x" << checkxi + 1 << "/" << "Xk" << '\n';

                zjcj.clear();
                pivoti = 0, pivotj = 0;
                minzjcj = 100100;
                minxbyr = 100100;

                for (int i = 0; i < c.size(); i++)
                {
                    double val = findzj(c, a, basicindex, i) - c[i];
                    zjcj.push_back(val);
                    if (minzjcj > val)
                    {
                        pivotj = i;
                        minzjcj = val;
                    }
                }

                xbyrvec.clear();

                basicsolution.clear();

                printMidLeft(basicindex, c, a, b, repeat, checkxi, pivotj, xbyrvec, minxbyr, pivoti, basicsolution);

                // for(auto i : basicsolution){
                //     cout<<"("<<i.first<<","<<i.second<<") ";
                // }cout<<endl;

                printLowerLine(zjcj);

                repeat = checkDegeracy(minxbyr, xbyrvec);

                // cout<<"\n\nrepeat \n";
                // for(int i=0 ; i<repeat.size() ; i++){
                //     cout<< repeat[i] <<"  ";
                // }cout<<"\n" ;
            }
        }

        //---->> if repeat off <<---

        cout << "\n\n";
        cout << "pivot ele (" << pivoti + 1 << " , " << pivotj + 1 << ") \n\n";

        if (minzjcj < 0)
        {
            basicindex.erase(basicindex.begin() + pivoti);
            basicindex.insert(basicindex.begin() + pivoti, pivotj);
        }

        // sort(basicindex.begin(),basicindex.end());
        // map<int, int> solmain;
        if (minzjcj >= 0 || abs(minzjcj) < 0.01)
        {

            if (!fhase1Complte)
            {
                return;
            }

            cout << GREEN << "ANS IS : \n"
                 << RESET;
            for (int i = 0; i < basicindex.size(); i++)
            {
                // cout << "X" << basicindex[i] + 1 << " = " << b[i] << endl;
                if (b[i] == 0)
                {
                    degeneracyType1 = true;
                }
                // solmain[basicindex[i]] = b[i];
                if (find(artificialval.begin(), artificialval.end(), basicindex[i]) != artificialval.end())
                {
                    cout << RED << "X" << basicindex[i] + 1 << " = " << b[i] << endl
                         << RESET;
                    infeasiblesol = true;
                    Nosol = true;
                }
                else
                {
                    cout << GREEN << "X" << basicindex[i] + 1 << " = " << b[i] << endl
                         << RESET;
                }

                double fractionalPart = b[i] - floor(b[i]);
                cout << fractionalPart << ' ';
            }

            int max_fractonal = 0;
            int max_frac_at = -1;
            int max_index_in_a = -1;

            for (int i = 0; i < basicindex.size(); i++)
            {
                if (basicindex[i] < numberOfVariabe)
                {
                    double fractionalPart = b[i] - floor(b[i]);
                    cout << fractionalPart << ' ';
                    if (max_fractonal < fractionalPart)
                    {
                        max_fractonal = fractionalPart;
                        max_frac_at = basicindex[i];
                        max_index_in_a = i;
                    }
                }
            }

            cout << "\nmax_fractonal " << max_fractonal << " max_frac_at " << max_frac_at << " max_index_in_a " << max_index_in_a << endl;

            if (max_frac_at == -1)
            {
                cout << GREEN << "\n ans is over ans \n"
                     << RESET;
            }
            else
            {
                cout << YELLOW << "\n\nwe need to use grmorian cutting variable \n\n"
                     << RESET;
                vector<vector<double>> new_a(a.size() + 1, vector<double>(a[0].size(), 0));
                for (int i = 0; i < a.size(); i++)
                {
                    for (int j = 0; j < a[0].size(); j++)
                    {
                        new_a[i][j] = a[i][j];
                    }
                }
                for (int j = 0; j < a[0].size(); j++)
                {
                    double fractionalPart = a[max_index_in_a][j] - floor(a[max_index_in_a][j]);
                    new_a[new_a.size()][j] = fractionalPart;
                }

                cout << "\nmax_fractonal " << max_fractonal << " max_frac_at " << max_frac_at << " max_index_in_a " << max_index_in_a << endl;

                for (int i = 0; i < new_a.size(); i++)
                {
                    for (int j = 0; j < new_a[0].size(); j++)
                    {
                        new_a[i][j] = a[i][j];
                    }
                }
            }

            if (Nosol)
            {
                cout << RED << "\n\nNon-existing feasible solution solution \n\n"
                     << RESET;
            }
            double ovalue = 0;

            for (int i = 0; i < numberOfVariabe; i++)
            {
                ovalue += basicsolution[i] * c[i];
                // cout << "\n" << i << ' ' << basicsolution[i] << ' ' << c[i] << "\n";
            }
            saveAlternatesol.push_back(basicsolution);
            if (infeasiblesol)
            {
                cout << RED << "optimal value : " << RESET;
                // if(!infeasiblesol)
                if (minTypeOptimalSol)
                {
                    cout << RED << -1 * ovalue << RESET << endl;
                }
                else
                {
                    cout << RED << ovalue << RESET << endl;
                }
            }
            else
            {
                cout << GREEN << "optimal value : " << RESET;
                // if(!infeasiblesol)
                if (minTypeOptimalSol)
                {
                    cout << GREEN << -1 * ovalue << RESET << endl;
                }
                else
                {
                    cout << GREEN << ovalue << RESET << endl;
                }
            }

            if (isinf(ovalue))
            {
                cout << RED << "\n\nunbounded solution\n\n"
                     << RESET;
                UnboundedSol = true;
                return;
            }

            vector<int> zeroIndex = checkAlternateSol(zjcj, basicindex);

            if (!fhase1Complte || zeroIndex.size() == basicindex.size())
            {
                // cout << "CHaeck1n\n\n";
                return;
            }
            else
            {
                // cout << "\n\nCHekc check non basic show 0\n\n";
                int nonbasic = checkNonbasicShowzero(zeroIndex, basicindex, a, pivoti);
                if (nonbasic == -1)
                {

                    return;
                }
                alternateSol = true;
                // cout << "\n\nthe existence of alternate solution\n\n";
                vector<int> tem = basicindex;
                sort(tem.begin(), tem.end());

                basicindex.erase(basicindex.begin() + pivoti);
                basicindex.insert(basicindex.begin() + pivoti, nonbasic);
                pivotj = nonbasic;
                // pivoti = tem[tem.size()-1] ;
                i = 18;
                isnosol = false;
                if (saveAlternatesol.size() >= 2)
                {
                    return;
                }
            }
        }

        operation(pivoti, pivotj, a, b);
    }

    if (isnosol)
    {
        cout << RED << "\n\n   --->>> NO SOLUTION <<<--- \n\n"
             << RESET;
        Nosol = true;
    }
}

int main()
{
    int n, m;

    // Input number of variables and constraints
    cout << VIOLET << "Enter the number of variables: " << RESET;
    cin >> n;
    numberOfVariabe = n;
    cout << VIOLET << "Enter the number of constraints: " << RESET;
    cin >> m;
    numberOfEquation = m;

    vector<vector<double>> A(m, vector<double>(n));
    vector<double> b(m);
    vector<double> c(n);

    string sing;
    vector<int> equal, greater, less;

    cout << VIOLET << "Enter the coefficients of the constraints (one line per constraint, separated by spaces ) :" << RESET << endl;
    // cout << "Just like Ai1 Ai2 .... Ain  bi \n";
    for (int i = 0; i < m; ++i)
    {
        cout << VIOLET << "Enter A(" << i + 1 << " , j) : \n"
             << RESET;
        for (int j = 0; j < n; ++j)
        {
            cin >> A[i][j];
        }

        cout << VIOLET << "Enter b(i) : " << RESET;
        cin >> b[i]; // Right-hand side of the constraint

        cout << VIOLET << "Enter sing between A(" << i + 1 << " , j)" << "and b(i) : " << RESET;
        cin >> sing;
        if (sing == "=" || sing == "==")
        {
            equal.push_back(i);
        }
        else if (sing == "<" || sing == "<=")
        {
            less.push_back(i);
        }
        else if (sing == ">" || sing == ">=")
        {
            greater.push_back(i);
        }
    }

    cout << VIOLET << "Enter the coefficients of the objective function (one line, separated by spaces):" << RESET << endl;
    for (int j = 0; j < n; ++j)
    {
        cin >> c[j];
    }
    cout << VIOLET << "Objective function is Min type(0) Or Max type(1) : " << RESET;
    int minmax;
    cin >> minmax;

    // cout<<"mention given Xi is resited or unresi"

    if (minmax == 0)
    {
        minTypeOptimalSol = true;
        for (int i = 0; i < c.size(); i++)
        {
            c[i] *= -1;
        }
    }

    // Convert to standard form
    vector<vector<double>> A_standard(m, vector<double>(n + equal.size() + 2 * (greater.size()) + less.size(), 0));
    vector<double> b_standard(m);
    vector<double> c_standard(n + equal.size() + 2 * (greater.size()) + less.size(), 0);

    int valfill = n;
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            A_standard[i][j] = A[i][j];
        }
        if (find(less.begin(), less.end(), i) != less.end())
        {
            A_standard[i][valfill] = 1; // Slack variable
            valfill++;
        }
        else if (find(greater.begin(), greater.end(), i) != greater.end())
        {
            A_standard[i][valfill] = -1;
            valfill++;
        }
        b_standard[i] = b[i];
    }

    for (int i = 0; i < m; i++)
    {
        if (find(greater.begin(), greater.end(), i) != greater.end())
        {
            A_standard[i][valfill] = 1;
            c_standard[valfill] = -1 * M;
            artificialval.push_back(valfill);
            valfill++;
        }
        if (find(equal.begin(), equal.end(), i) != equal.end())
        {
            A_standard[i][valfill] = 1;
            c_standard[valfill] = -1 * M;
            artificialval.push_back(valfill);
            valfill++;
        }
    }

    for (int j = 0; j < n; ++j)
    {
        c_standard[j] = c[j];
    }

    vector<int> basic_index;

    // for (int i = c_standard.size() - 1; i >= c_standard.size() - m; i--)
    // {
    //     basic_index.push_back(i);
    // }

    for (int i = 0; i < A_standard.size(); i++)
    {
        int j = A_standard[i].size() - 1;
        while (j >= 0 && A_standard[i][j] == 0)
        {
            j--;
        }
        basic_index.push_back(j);
    }

    // reverse(basic_index.begin(), basic_index.end());
    // basic_index = {2, 3, 5};
    // auto tem1 = basic_index.back() ;
    // basic_index.pop_back();
    // basic_index.insert(basic_index.begin() , tem1) ;
    // tem1 = basic_index.back() ;
    // basic_index.pop_back();
    // basic_index.insert(basic_index.begin() , tem1) ;
    cout << "\n";
    vector<double> cprime(c_standard.size(), 0);
    for (int i = 0; i < c_standard.size(); i++)
    {
        if (c_standard[i] == -1 * M)
        {
            cprime[i] = -1;
        }
    }
    if (artificialval.size() > 0)
    {
        printtable(A_standard, b_standard, cprime, basic_index);
        for (auto bi : basic_index)
        {
            if (find(artificialval.begin(), artificialval.end(), bi) != artificialval.end())
            {
                cout << RED << "\nNon-Existing feasible Solution \n\n"
                     << RESET;
                return 0;
            }
        }
        int needToRemove = artificialval.size();

        for (int i = 0; i < needToRemove; i++)
        {
            c_standard.pop_back();
            for (int j = 0; j < A_standard.size(); j++)
            {
                A_standard[j].pop_back();
            }
        }
    }
    fhase1Complte = true;
    printtable(A_standard, b_standard, c_standard, basic_index);

    if (degeneracy || degeneracyType1)
    {
        cout << "\n\nSolution show degeneracy \n\n";
    }
    if (alternateSol)
    {
        cout << "\n\nthe existence of alternate solution\n\n";

        cout << "Variables   First Sol.   Second Sol.   General Solution\n";
        for (int i = 0; i < c_standard.size(); i++)
        {
            cout << "   X" << i + 1 << "        " << saveAlternatesol[0][i] << "           " << saveAlternatesol[1][i] << "        " << "   x" << i + 1 << " = (" << saveAlternatesol[0][i] << ")*$ + (" << saveAlternatesol[1][i] << ")*( 1 - $ )\n";
        }
    }
    if (infeasiblesol && !UnboundedSol)
    {
        cout << RED << "\nNon-Existing feasible Solution \n\n"
             << RESET;
    }

    return 0;
}

// ass 6 page 8 example 23 doubte