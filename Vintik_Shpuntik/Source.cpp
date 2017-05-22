#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <iomanip> 
#define det(a,b,c,d)  (a*d-b*c)
const double EPS = 1E-9;
using namespace std;
int N,num=0,vin_r,shpun_r;
//vector < vector < pair<int, double> > > graph_main(1000);
vector <pair<pair<int, int>, pair<int, int>>> data_for_graph;
//vector <pair<pair<int, int>, pair<int, int>>> vec_otr;
vector <pair<pair<int, int>, pair<int, int>>> vec_help;
vector < vector < pair<int, double> > > graph_main(1000);
pair<int, int> vint, shpun;

struct pt {
	double x, y;

	bool operator< (const pt & p) const {
		return x < p.x - EPS || abs(x - p.x) < EPS && y < p.y - EPS;
	}
};

struct line {
	double a, b, c;

	line() {}
	line(pt p, pt q) {
		a = p.y - q.y;
		b = q.x - p.x;
		c = -a * p.x - b * p.y;
		norm();
	}

	void norm() {
		double z = sqrt(a*a + b*b);
		if (abs(z) > EPS)
			a /= z, b /= z, c /= z;
	}

	double dist(pt p) const {
		return a * p.x + b * p.y + c;
	}
};

inline bool intersect_1d(double a, double b, double c, double d) {
	if (a > b)  swap(a, b);
	if (c > d)  swap(c, d);
	return max(a, c) <= min(b, d) + EPS;
}

inline bool betw(double l, double r, double x) {
	return min(l, r) <= x + EPS && x <= max(l, r) + EPS;
}

pair<int,int> Cross_Line(int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4) {
	if (!intersect_1d(x1, x2, x3, x4) || !intersect_1d(y1, y2, y3, y4))
		return make_pair(12345,1235);
	pt a,b,c,d,left,right;
	a.x = x1;
	a.y = y1;
	b.x = x2;
	b.y = y2;
	c.x = x3;
	c.y = y3;
	d.x = x4;
	d.y = y4;
	line m(a, b);
	line n(c, d);
	double zn = det(m.a, m.b, n.a, n.b);
	if (abs(zn) < EPS) {
		if (abs(m.dist(c)) > EPS || abs(n.dist(a)) > EPS)
			return make_pair(12345,12345);
		if (b < a)  swap(a, b);
		if (d < c)  swap(c, d);
		left = max(a, c);
		right = min(b, d);
		return make_pair(left.x,left.y);
	}
	else {
		left.x = right.x = -det(m.c, m.b, n.c, n.b) / zn;
		left.y = right.y = -det(m.a, m.c, n.a, n.c) / zn;
		if(betw(a.x, b.x, left.x) && betw(a.y, b.y, left.y) && betw(c.x, d.x, left.x) && betw(c.y, d.y, left.y))
			return make_pair(left.x, left.y);
		else return make_pair(12345, 12345);
	}
}

pair<int, int> Cross_Line_main(int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4)
{
	pair<int, int> p1 = Cross_Line(x1, y1, x2, y2, x3, y3, x4, y4);
	pair<int, int> p2 = Cross_Line(x3, y3, x4, y4, x1, y1, x2, y2);
	if (p1.first != 12345 && p1.second != 12345)
		return p1;
	else if (p2.first != 12345 && p2.second != 12345)
		return p2;
	else return make_pair(12345, 12345);
}

double calculate_angle(const pair<int, int> &p1, const pair<int, int> &p2, const pair<int, int> &p3) //todo возможно не будет работать с отрицательными числами
{
	//if ((p1.first == p3.first && p1.second != p3.second) || (p1.first != p3.first && p1.second == p3.second))
	if ((p1.first==p2.first && p1.first==p3.first) || (p1.second == p2.second && p1.second == p3.second))
		return 0;
	pair<int, int> vec1(p1.first - p2.first, p1.second - p2.second);
	pair<int, int> vec2(p3.first - p2.first, p3.second - p2.second);
	const double PI = 3.14159265;
	const double length_vec1 = sqrt((vec1.first*vec1.first) + (vec1.second*vec1.second));
	const double length_vec2 = sqrt((vec2.first*vec2.first) + (vec2.second*vec2.second));
	const double dotProduct = vec1.first*vec2.second + vec1.second*vec2.first;
	const double cosTheta = dotProduct / (length_vec1*length_vec2);
	const double res = (acos(cosTheta)*180.0) / PI;
	if (res == 0)
		return 90.0;
	else
		return res;
}

void read_from_file()
{
	ifstream ifs("input.txt");
	ifs >> N;

	int x1,y1,x2,y2;

	for (int i = 0;i<N;++i)
	{
		ifs >> x1 >> y1 >> x2 >> y2;
		data_for_graph.emplace_back(make_pair(x1, y1), make_pair(x2, y2));
	}
	ifs >> x1 >> y1 >> x2 >> y2;
	vint = make_pair(x1, y1);
	shpun = make_pair(x2, y2);
	ifs.close();
}

void build_graph()
{
	bool flag = 0;
	vector < pair<int, double> > vec_h;
	
	//int num = 0;

	for (int i = 0;i < N - 1;++i)
	{
		int x11 = data_for_graph[i].first.first;
		int y11 = data_for_graph[i].first.second;
		int x12 = data_for_graph[i].second.first;
		int y12 = data_for_graph[i].second.second;
		for (int j = i + 1;j < N;++j)
		{
			int x21 = data_for_graph[j].first.first;
			int y21 = data_for_graph[j].first.second;
			int x22 = data_for_graph[j].second.first;
			int y22 = data_for_graph[j].second.second;
			pair<int, int> p = Cross_Line_main(x11, y11, x12, y12, x21, y21, x22, y22);
			if (p.first != 12345 && p.second != 12345)
			{
				if (x11 == p.first && y11 == p.second || x12 == p.first && y12 == p.second)
				{
					pair<pair<int, int>, pair<int, int>> p_h = make_pair(make_pair(x11, y11), make_pair(x12, y12));
					vec_help.emplace_back(p_h);
					num = num + 1;
				}
				else
					{
						pair<pair<int, int>, pair<int, int>> p_h1 = make_pair(make_pair(x11, y11), make_pair(p.first, p.second));
						pair<pair<int, int>, pair<int, int>> p_h2 = make_pair(make_pair(p.first, p.second), make_pair(x12, y12));
						vec_help.emplace_back(p_h1);
						vec_help.emplace_back(p_h2);
						num = num + 1;
						vec_help.erase(remove(vec_help.begin(), vec_help.end(), data_for_graph[i]), vec_help.end());
					}

				if (x21 == p.first && y21 == p.second || x22 == p.first && y22 == p.second)
				{
					pair<pair<int, int>, pair<int, int>> p_h = make_pair(make_pair(x21, y21), make_pair(x22, y22));
					vec_help.emplace_back(p_h);
					num = num + 1;
				}
				else
				{
					pair<pair<int, int>, pair<int, int>> p_h1 = make_pair(make_pair(x21, y21), make_pair(p.first, p.second));
					pair<pair<int, int>, pair<int, int>> p_h2 = make_pair(make_pair(p.first, p.second), make_pair(x22, y22));
					vec_help.emplace_back(p_h1);
					vec_help.emplace_back(p_h2);
					num = num + 1;
					vec_help.erase(remove(vec_help.begin(), vec_help.end(), data_for_graph[j]), vec_help.end());
				}
			}
		}
	}
	int x11 = data_for_graph[0].first.first;
	int y11 = data_for_graph[0].first.second;
	int x12 = data_for_graph[0].second.first;
	int y12 = data_for_graph[0].second.second;
	set<pair<pair<int, int>, pair<int, int>>> b(vec_help.begin(), vec_help.end());
	vector <pair<pair<int, int>, pair<int, int>>> vec_otr(b.begin(), b.end());
	if (N == 1)
	{
		num = 1;
		vec_otr.emplace_back(make_pair(x11, y11), make_pair(x12, y12));
	}
	//cout << num << endl;
	
	num = vec_otr.size();
	//cout << num << endl;
	double weight;
	//graph_main.assign();
	for (int i = 0;i < num;i++)
	{
		int x11 = vec_otr[i].first.first;
		int y11 = vec_otr[i].first.second;
		int x12 = vec_otr[i].second.first;
		int y12 = vec_otr[i].second.second;
		if (((vint.first >= x11 && vint.first <= x12) || (vint.first <= x11 && vint.first >= x12)) && ((vint.second >= y11 && vint.second <= y12) || (vint.second <= y11 && vint.second >= y12)))
			vin_r = i;
		if (((shpun.first >= x11 && shpun.first <= x12) || (shpun.first <= x11 && shpun.first >= x12)) && ((shpun.second >= y11 && shpun.second <= y12) || (shpun.second <= y11 && shpun.second >= y12)))
			shpun_r = i;
		for (int j = 0;j < num;j++)
		{
			if (i != j)
			{
				int x21 = vec_otr[j].first.first;
				int y21 = vec_otr[j].first.second;
				int x22 = vec_otr[j].second.first;
				int y22 = vec_otr[j].second.second;
				pair<int, int> p3 = Cross_Line_main(x11, y11, x12, y12, x21, y21, x22, y22);
				if (p3.first != 12345 && p3.second != 12345)
				{
					pair<int, int> p1, p2;
					p1 = make_pair(x11, y11);
					p2 = make_pair(x21, y21);

					if (p1 == p3)
						p1 = make_pair(x12, y12);
					if (p2 == p3)
						p2 = make_pair(x22, y22);
					weight = calculate_angle(p1, p3, p2);
					//vec_h.emplace_back(j, weight);
					graph_main[i].emplace_back(j, weight);
				}
			}
			
		}
		vec_h.clear();
	}



	int t = 0;
}

const int INF = 1000000000;

double deikstra(vector < vector < pair<int, double> > > g) {
	int n = num;
	int s = vin_r;
	vector<double> d(n, INF);
	vector<int>	p(n);

	d[s] = 0;
	vector<char> u(n);
	for (int i = 0; i<n; ++i) {
		int v = -1;
		for (int j = 0; j<n; ++j)
			if (!u[j] && (v == -1 || d[j] < d[v]))
				v = j;
		if (d[v] == INF)
			break;
		u[v] = true;

		for (size_t j = 0; j<g[v].size(); ++j) {
			int to = g[v][j].first;
			double len = g[v][j].second;
			if (d[v] + len < d[to]) {
				d[to] = d[v] + len;
				p[to] = v;
			}
		}
	}
	return d[shpun_r];
}

int main()
{
	read_from_file();
	//cout << N << endl;
	//cout << vint.second << endl;
	//find_graph(vec_2_pairs_hor, vec_2_pairs_ver);
	//pair<int,int> res= Cross_Line_main(0, 1, 1, 1, 0, 0, -1, 1);
	double res;
	build_graph();
	if (num == 0)
		res = -1;
	else
		res = deikstra(graph_main);
	//cout << res << endl;
	//cout << CrossLine(0, 0, 1, 0, 0, 0, 0, 2).first << " " << CrossLine(0, 0, 1, 0, 0, 0, 0, 2).second << endl;
//	cout << res.first << " " << res.second << endl;
	//int res = 0;
	//res = 3.1617384;
	//cout << setprecision(5) << res << endl;
	ofstream ofs("output.txt");
		ofs << setprecision(5) << res << endl;
	ofs.close();


	return 0;
}