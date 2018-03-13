/**
 * @file Profiler.h
 * @author Fabjan
 * @date 9/02/2015
 * @version 1.0
 * 
 * Customizable profiler
 */

#pragma once

#include <emotion/3rdparty/Timer.h>

#include <map>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

#define PRUN(X)\
g_prof.run(X)

#define PSTOP(X)\
g_prof.stop(X)

struct ProfilingElement {

	ProfilingElement() {
		count = 0;
		avg = 0.0;
		max = 0;
		min = 1e9;
	}

	void update() {
		last = t.getTime().milli;
		avg = (avg * count + last) / (++count);
		if (last > max) max = last;
		if (last < min) min = last;
	}

	void start() {
		t.start();
	}

	Timer t;
	int count;
	double avg;
	double last;
	double max;
	double min;
};

class Profiler {
public:

	Profiler() {
	};

	void run(const std::string& name) {
		auto it = m_prof.find(name);
		if (it != m_prof.end()) {
			it->second.start();
		} else {
			m_prof[name] = *(new ProfilingElement());
			m_prof[name].t.start();
		}
	}

	void stop(const std::string& name) {
		auto it = m_prof.find(name);
		if (it != m_prof.end()) {
			it->second.update();
		}
	}

	double get(const std::string& name) {
		return m_prof[name].last;
	}

	double reset(const std::string& name) {
		m_prof[name].start();
	}

	void printInfo(double scale = -1) {

		std::vector<std::pair<std::string, ProfilingElement> > pairs(m_prof.begin(), m_prof.end());

		std::sort(pairs.begin(), pairs.end(), [ = ](const std::pair<std::string, ProfilingElement>& a, const std::pair<std::string, ProfilingElement>& b){
			return a.second.avg > b.second.avg;
		}
		);

		int c = 0;
		for (auto it = pairs.begin(); it != pairs.end(); ++it) {
			std::cout << c++ << ": " << it->first << "\t>\tavg:\t" << it->second.avg << "\tcount:\t" << it->second.count;
			std::cout << "\tmin:\t" << it->second.min << "\tmax:\t" << it->second.max << "\tlast:\t" << it->second.last;
			if (scale > 0) {
				std::cout << "\t%:\t" << 100 * it->second.avg / scale;
			}
			std::cout << std::endl;
		}
	}
private:
	std::map<std::string, ProfilingElement> m_prof;
};

extern Profiler g_prof;
