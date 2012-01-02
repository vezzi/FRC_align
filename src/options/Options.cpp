/*
 *  This file is part of rNA.
 *  Copyright (c) 2011 by Cristian Del Fabbro <delfabbro@appliedgenomics.org>,
 *  Francesco Vezzi <vezzi@appliedgenomics.org>,
 *  Alexandru Tomescu <alexandru.tomescu@uniud.it>, and
 *  Alberto Policriti <policriti@uniud.it>
 *
 *   rNA is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.

 *   rNA is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.

 *   You should have received a copy of the GNU General Public License
 *   along with rNA.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "Options.h"

namespace options {

void Options::set_defaults() {
	program_mode = program_unknown;

	auto_errors = false;
	errors_rate = 15;
	common_errors_allowed = 0;

	force_fastqformat = false;

	contamination_check = false;

	quality_check = true;
	trim = true;
	min_phred_value_CLC = 20;
	min_mean_quality = 20;
	min_size = 25;

	seed_sizes = 10;
	seed_errors = 1;
	max_gap = 10;
	gap = false;

	sample = string("no_sample_specified");
	bam_format = false;

	threads_number = 1;

	paired_ends = false;

	insert_size_check = false;

	k = 12;
	blockLength = 10;

	delta = 0;
	indels = false;
	indels_max_value = 5;

	// output options
	printAll = false;
	toBePrinted = UINT_MAX;

	gui_output = false;

}

}
