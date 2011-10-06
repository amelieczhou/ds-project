// TeleNetwork_simulation.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <queue>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <math.h>
#include <sstream>
//
using namespace std;


int Cell[20];
//int Cell9[20];
int BlockedCalls = 0;
int droppedCalls = 0;
int TotalCalls = 0;

class Event
{
public:
	double time;
	virtual void handler() = 0;//pure virtual function
	virtual string print_name() = 0;
	//virtual void handler() {cout<<"base class";}
};
class initiater:public Event //program only those differ from or extend the base class
{
public:
	double velocity;
	int station;
	double duration;
	double position;

	initiater(double t, double v, int s, double d, double p){time = t; velocity = v; station = s; duration = d; position = p; }
	string print_name() {stringstream ss; ss<<"i'm initiater";return ss.str();}
	virtual void handler();
};
class terminator:public Event
{
public:
	int station;

	terminator(double t, int s) {time = t; station = s; }
	string print_name(){stringstream ss;ss<<"i'm terminator";return ss.str();}
	virtual void handler();
};
class handover:public Event
{
public:
	double velocity;
	int station;
	double duration;

	handover(double t, double v, int s, double d) {time = t; velocity = v; station = s; duration = d; }
	string print_name(){stringstream ss; ss<<"i'm handover"; return ss.str();}
	virtual void handler();
};

struct cmp{
    bool operator() (Event* a, Event* b ){
		double t1 = a->time;
		double t2 = b->time;
		return t1 > t2; }
};
/*struct cmp{
	bool operator() (Event a, Event b){
		return a.time>b.time;
	}
}*/
priority_queue<Event*, vector<Event*>, cmp> sortedQueue;
//priority_queue<Event, vector<Event>, cmp> sortedQueue;
void initiater::handler()
{
	if(Cell[station] > 0)
	{
		Cell[station] -= 1;
		//position = station * 2;
		double endPosition = position + velocity * duration / 3600;//km per h
		int endStation = endPosition / 2;
		double travelledDuration = ((station + 1)*2 - position) / velocity * 3600;//sec
		if ( station == endStation)
		{
			terminator* ter_event = new terminator(time + duration, endStation);
			sortedQueue.push(ter_event);
		}			
		else if (endStation > 19)
		{
			terminator* ter_event = new terminator(time + travelledDuration, station);
			sortedQueue.push(ter_event);
		}
		else
		{
			handover* hand_event = new handover(time + travelledDuration, velocity, station + 1, duration - travelledDuration);
			sortedQueue.push(hand_event);
		}
	}
	else
		BlockedCalls += 1;
	TotalCalls += 1;
}
void terminator::handler()
{
	Cell[station] += 1;
}
void handover::handler()
{
	Cell[station - 1] +=1;
	
	if(Cell[station] > 0)
	{
		Cell[station] -= 1;
		double position = station * 2;
		double endPosition = position + velocity * duration / 3600;
		int endStation = endPosition / 2;
		double travelledDuration = ((station + 1)*2 - position) / velocity * 3600;
		if ( station == endStation)
		{
			terminator* ter_event = new terminator(time + duration, endStation);
			sortedQueue.push(ter_event);
		}			
		else if (endStation > 19)
		{
			terminator* ter_event = new terminator(time + travelledDuration, station);
			sortedQueue.push(ter_event);
		}
		else
		{
			handover* hand_event = new handover(time + travelledDuration, velocity, station + 1, duration - travelledDuration);
			sortedQueue.push(hand_event);
		}
	}
	else
		droppedCalls += 1;
}
int _tmain(int argc, _TCHAR* argv[])
{
	ofstream ofs_result;
	ofs_result.open("blocked&dropped rate under statistic case.txt",ios::out | ios::app);
	
	//generate random numbers with diff distribution
	/*srand( (unsigned)time( NULL ) );
	int n =500 + iter*500; //num of random numbers
	vector<double> time(n,0);
	vector<double> position(n,0);
	vector<double> duration(n,0);
	vector<double> velocity(n,0);
	for(int i=0; i<n; i++)
	{
		//uniform
		time[i] = (double)rand() / (RAND_MAX + 1);
		//triangular
		if(time[i] <= 0.75) position[i] = pow(1200*time[i],0.5);
		else position[i] = 40 - 20*pow(1-time[i], 0.5);
		//exponential
		duration[i] = -160 * log(time[i]);

		if(i != 0)	time[i] = time[i - 1] + time[i];

	}
	//normal distribution
	double _n = 12;//recommended
	for(int i=0; i<n; i++)
	{	double sum=0;
		for(int i=0; i<_n; i++)
			sum += (double)rand() / (RAND_MAX + 1);
		velocity[i] = 100 + 9*((sum - _n/2)/pow(_n/12, 0.5));
	}

	ofstream uni_ofs;
	uni_ofs.open("uniform distribution.txt",ios::out);
	ofstream exp_ofs;
	exp_ofs.open("exponential distribution.txt",ios::out);
	ofstream tri_ofs;
	tri_ofs.open("triangular distribution.txt",ios::out);
	ofstream nor_ofs;
	nor_ofs.open("normal distribution.txt",ios::out);
	for(int i=0; i<n; i++) 
	{
		uni_ofs<<time[i]<<endl;
		exp_ofs<<duration[i]<<endl;
		tri_ofs<<position[i]<<endl;
		nor_ofs<<velocity[i]<<endl;
	}*/

	int sample = 10000;
	//event handling process
	for(int i=0; i<20; i++) 
	{
		Cell[i] = 10;
		//Cell9[i] = 1;
	}
	
	//pending event list
	ifstream myfile("PCS_TEST_DETERMINSTIC.txt");
	if(myfile.is_open())
	{
		string line;
		getline(myfile,line);
		for(int i=0; i<sample; i++)
		{
			double t,v,d,p;
			int s;
			getline(myfile,line);
			char pattern = ',';
			int pos1 = line.find(pattern, 0);
			int pos2 = line.find(pattern, pos1 + 1);
			int pos3 = line.find(pattern, pos2 + 1);
			int pos4 = line.find(pattern, pos3 + 1);
			string s1 = line.substr(pos1 + 1,pos2 - pos1 - 1);
			t = stod(s1);
			string s2 = line.substr(pos2 + 1, pos3 - pos2 - 1);
			s = stoi(s2) - 1;//start from 0
			string s3 = line.substr(pos3 + 1, pos4 - pos3 - 1);
			d = stod(s3);
			string s4 = line.substr(pos4 + 1, line.length() - pos4 - 1);
			v = stod(s4);
			p = s * 2 + 1;//initialized at the centre of base stations
			bool c = 1; //channelId of initiallization is 1
			initiater* callEvent = new initiater(t,v,s,d,p);
			sortedQueue.push(callEvent);
		}
	}
	ofstream ofs_event("temp for xiaosong.txt");

	clock_t curr_t, end_t;
	curr_t = clock();
	while(!sortedQueue.empty())
	{
		Event* currEvent = sortedQueue.top();
		//Event currEvent = sortedQueue.top();
		currEvent->handler();
		ofs_event<<currEvent->print_name()<<"\t";
		ofs_event<<currEvent->time<<endl;			

		//cout<<currEvent->time<<endl;
		//currEvent.handler();
		sortedQueue.pop();
	}
	end_t = clock();
	double run_duration = (double)(end_t - curr_t)/CLOCKS_PER_SEC;//25.484s under debug mode
	double blocked_rate = (double)BlockedCalls / TotalCalls;//scheme2, 10000, 57; s1, 10000, 18
	double dropped_rate = (double)droppedCalls / TotalCalls;//scheme2, 10000, 11; s1, 10000, 18

	cout<<"run_duration: "<<run_duration<<", blocked_rate: "<<blocked_rate<<", dropped_rate: "<<dropped_rate<<endl;
	
	ofs_result<<blocked_rate<<'\t'<<dropped_rate<<endl;
	
	return 0;
}

