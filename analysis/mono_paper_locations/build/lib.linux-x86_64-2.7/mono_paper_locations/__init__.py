# before running analysis, set the environment variable co2_paper_locations to something that corresponds to this document
# e.g. export mono_paper_locations=analysiscavetech-original

#import data_locations
import figure_template_locations

import os
mono_paper_locations_profile = os.environ['mono_paper_locations']
print 'USING PROFILE: ', mono_paper_locations_profile

if mono_paper_locations_profile == 'analysiscavetech-organized':
	import data_locations_analysiscavetech_organized as data_locations
elif mono_paper_locations_profile == 'analysiscavetech-portable':
	import data_locations_analysiscavetech_portable as data_locations
