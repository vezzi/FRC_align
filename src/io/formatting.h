/*
 * formatting.h
 *
 *  Created on: 8-gen-2009
 *      Author: cdf
 */

#ifndef FORMATTING_H_
#define FORMATTING_H_

#include <sstream>
#include <iostream>

static inline std::ostream & print_formatted_time (std::ostream & os, double time) {
	double t = time / (3600*24);
	if (t > 1) {
		os << floor(t) << 'g';
		time -= floor(t)*(3600*24);
	}
	t = time / 3600;
	if (t > 1) {
		os << floor(t) << 'h';
		time -= floor(t)*3600;
	}
	t = time / 60;
	if (t > 1) {
		os << floor(t) << 'm';
		time -= floor(t)*60;
	}
	os << floor(time) << 's';
	return os;
}

#define TRUNC(x,d) ( (d <= -1) ? x : round(x*d)/d )

static inline std::string formatted_unit (double value, int digits = 1) {
	std::stringstream ss;
	double v = value / (1000.0*1000*1000*1000);
	long d;
	if (digits < 0)
		d = -1;
	else
		d = 10^digits;
	if (v >= 1) {
		ss << TRUNC(v,d) << 'T';
		return ss.str();
	}

	v = value / (1000.0*1000*1000);
	if (v >= 1) {
		ss << TRUNC(v,d) << 'G';
		return ss.str();
	}

	v = value / (1000.0*1000);
	if (v >= 1) {
		ss << TRUNC(v,d) << 'M';
		return ss.str();
	}

	v = value / (1000.0);
	if (v >= 1) {
		ss << TRUNC(v,d) << 'K';
		return ss.str();
	}

	ss << TRUNC(value,d);
	return ss.str();
}

#endif /* FORMATTING_H_ */
