//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Copyright (C) 2009  Sebastian Mach (*1983)
// * phresnel/at/gmail/dot/com
// * http://phresnel.org
// * http://picogen.org
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifndef RENDERER_HH_INCLUDED_20090503
#define RENDERER_HH_INCLUDED_20090503



#include <cmath>

template <int size> struct MakeIndex {
        size_t operator () (int x, int y) const {
                return x + y * size;
        }
};




namespace Morton {
        static const unsigned int MortonTable256[] = 
        {
          0x0000, 0x0001, 0x0004, 0x0005, 0x0010, 0x0011, 0x0014, 0x0015, 
          0x0040, 0x0041, 0x0044, 0x0045, 0x0050, 0x0051, 0x0054, 0x0055, 
          0x0100, 0x0101, 0x0104, 0x0105, 0x0110, 0x0111, 0x0114, 0x0115, 
          0x0140, 0x0141, 0x0144, 0x0145, 0x0150, 0x0151, 0x0154, 0x0155, 
          0x0400, 0x0401, 0x0404, 0x0405, 0x0410, 0x0411, 0x0414, 0x0415, 
          0x0440, 0x0441, 0x0444, 0x0445, 0x0450, 0x0451, 0x0454, 0x0455, 
          0x0500, 0x0501, 0x0504, 0x0505, 0x0510, 0x0511, 0x0514, 0x0515, 
          0x0540, 0x0541, 0x0544, 0x0545, 0x0550, 0x0551, 0x0554, 0x0555, 
          0x1000, 0x1001, 0x1004, 0x1005, 0x1010, 0x1011, 0x1014, 0x1015, 
          0x1040, 0x1041, 0x1044, 0x1045, 0x1050, 0x1051, 0x1054, 0x1055, 
          0x1100, 0x1101, 0x1104, 0x1105, 0x1110, 0x1111, 0x1114, 0x1115, 
          0x1140, 0x1141, 0x1144, 0x1145, 0x1150, 0x1151, 0x1154, 0x1155, 
          0x1400, 0x1401, 0x1404, 0x1405, 0x1410, 0x1411, 0x1414, 0x1415, 
          0x1440, 0x1441, 0x1444, 0x1445, 0x1450, 0x1451, 0x1454, 0x1455, 
          0x1500, 0x1501, 0x1504, 0x1505, 0x1510, 0x1511, 0x1514, 0x1515, 
          0x1540, 0x1541, 0x1544, 0x1545, 0x1550, 0x1551, 0x1554, 0x1555, 
          0x4000, 0x4001, 0x4004, 0x4005, 0x4010, 0x4011, 0x4014, 0x4015, 
          0x4040, 0x4041, 0x4044, 0x4045, 0x4050, 0x4051, 0x4054, 0x4055, 
          0x4100, 0x4101, 0x4104, 0x4105, 0x4110, 0x4111, 0x4114, 0x4115, 
          0x4140, 0x4141, 0x4144, 0x4145, 0x4150, 0x4151, 0x4154, 0x4155, 
          0x4400, 0x4401, 0x4404, 0x4405, 0x4410, 0x4411, 0x4414, 0x4415, 
          0x4440, 0x4441, 0x4444, 0x4445, 0x4450, 0x4451, 0x4454, 0x4455, 
          0x4500, 0x4501, 0x4504, 0x4505, 0x4510, 0x4511, 0x4514, 0x4515, 
          0x4540, 0x4541, 0x4544, 0x4545, 0x4550, 0x4551, 0x4554, 0x4555, 
          0x5000, 0x5001, 0x5004, 0x5005, 0x5010, 0x5011, 0x5014, 0x5015, 
          0x5040, 0x5041, 0x5044, 0x5045, 0x5050, 0x5051, 0x5054, 0x5055, 
          0x5100, 0x5101, 0x5104, 0x5105, 0x5110, 0x5111, 0x5114, 0x5115, 
          0x5140, 0x5141, 0x5144, 0x5145, 0x5150, 0x5151, 0x5154, 0x5155, 
          0x5400, 0x5401, 0x5404, 0x5405, 0x5410, 0x5411, 0x5414, 0x5415, 
          0x5440, 0x5441, 0x5444, 0x5445, 0x5450, 0x5451, 0x5454, 0x5455, 
          0x5500, 0x5501, 0x5504, 0x5505, 0x5510, 0x5511, 0x5514, 0x5515, 
          0x5540, 0x5541, 0x5544, 0x5545, 0x5550, 0x5551, 0x5554, 0x5555
        };
        static const unsigned int MortonTable256shl1[] = 
        {
          0x0000 << 1, 0x0001 << 1, 0x0004 << 1, 0x0005 << 1, 0x0010 << 1, 0x0011 << 1, 0x0014 << 1, 0x0015 << 1, 
          0x0040 << 1, 0x0041 << 1, 0x0044 << 1, 0x0045 << 1, 0x0050 << 1, 0x0051 << 1, 0x0054 << 1, 0x0055 << 1, 
          0x0100 << 1, 0x0101 << 1, 0x0104 << 1, 0x0105 << 1, 0x0110 << 1, 0x0111 << 1, 0x0114 << 1, 0x0115 << 1, 
          0x0140 << 1, 0x0141 << 1, 0x0144 << 1, 0x0145 << 1, 0x0150 << 1, 0x0151 << 1, 0x0154 << 1, 0x0155 << 1, 
          0x0400 << 1, 0x0401 << 1, 0x0404 << 1, 0x0405 << 1, 0x0410 << 1, 0x0411 << 1, 0x0414 << 1, 0x0415 << 1, 
          0x0440 << 1, 0x0441 << 1, 0x0444 << 1, 0x0445 << 1, 0x0450 << 1, 0x0451 << 1, 0x0454 << 1, 0x0455 << 1, 
          0x0500 << 1, 0x0501 << 1, 0x0504 << 1, 0x0505 << 1, 0x0510 << 1, 0x0511 << 1, 0x0514 << 1, 0x0515 << 1, 
          0x0540 << 1, 0x0541 << 1, 0x0544 << 1, 0x0545 << 1, 0x0550 << 1, 0x0551 << 1, 0x0554 << 1, 0x0555 << 1, 
          0x1000 << 1, 0x1001 << 1, 0x1004 << 1, 0x1005 << 1, 0x1010 << 1, 0x1011 << 1, 0x1014 << 1, 0x1015 << 1, 
          0x1040 << 1, 0x1041 << 1, 0x1044 << 1, 0x1045 << 1, 0x1050 << 1, 0x1051 << 1, 0x1054 << 1, 0x1055 << 1, 
          0x1100 << 1, 0x1101 << 1, 0x1104 << 1, 0x1105 << 1, 0x1110 << 1, 0x1111 << 1, 0x1114 << 1, 0x1115 << 1, 
          0x1140 << 1, 0x1141 << 1, 0x1144 << 1, 0x1145 << 1, 0x1150 << 1, 0x1151 << 1, 0x1154 << 1, 0x1155 << 1, 
          0x1400 << 1, 0x1401 << 1, 0x1404 << 1, 0x1405 << 1, 0x1410 << 1, 0x1411 << 1, 0x1414 << 1, 0x1415 << 1, 
          0x1440 << 1, 0x1441 << 1, 0x1444 << 1, 0x1445 << 1, 0x1450 << 1, 0x1451 << 1, 0x1454 << 1, 0x1455 << 1, 
          0x1500 << 1, 0x1501 << 1, 0x1504 << 1, 0x1505 << 1, 0x1510 << 1, 0x1511 << 1, 0x1514 << 1, 0x1515 << 1, 
          0x1540 << 1, 0x1541 << 1, 0x1544 << 1, 0x1545 << 1, 0x1550 << 1, 0x1551 << 1, 0x1554 << 1, 0x1555 << 1, 
          0x4000 << 1, 0x4001 << 1, 0x4004 << 1, 0x4005 << 1, 0x4010 << 1, 0x4011 << 1, 0x4014 << 1, 0x4015 << 1, 
          0x4040 << 1, 0x4041 << 1, 0x4044 << 1, 0x4045 << 1, 0x4050 << 1, 0x4051 << 1, 0x4054 << 1, 0x4055 << 1, 
          0x4100 << 1, 0x4101 << 1, 0x4104 << 1, 0x4105 << 1, 0x4110 << 1, 0x4111 << 1, 0x4114 << 1, 0x4115 << 1, 
          0x4140 << 1, 0x4141 << 1, 0x4144 << 1, 0x4145 << 1, 0x4150 << 1, 0x4151 << 1, 0x4154 << 1, 0x4155 << 1, 
          0x4400 << 1, 0x4401 << 1, 0x4404 << 1, 0x4405 << 1, 0x4410 << 1, 0x4411 << 1, 0x4414 << 1, 0x4415 << 1, 
          0x4440 << 1, 0x4441 << 1, 0x4444 << 1, 0x4445 << 1, 0x4450 << 1, 0x4451 << 1, 0x4454 << 1, 0x4455 << 1, 
          0x4500 << 1, 0x4501 << 1, 0x4504 << 1, 0x4505 << 1, 0x4510 << 1, 0x4511 << 1, 0x4514 << 1, 0x4515 << 1, 
          0x4540 << 1, 0x4541 << 1, 0x4544 << 1, 0x4545 << 1, 0x4550 << 1, 0x4551 << 1, 0x4554 << 1, 0x4555 << 1, 
          0x5000 << 1, 0x5001 << 1, 0x5004 << 1, 0x5005 << 1, 0x5010 << 1, 0x5011 << 1, 0x5014 << 1, 0x5015 << 1, 
          0x5040 << 1, 0x5041 << 1, 0x5044 << 1, 0x5045 << 1, 0x5050 << 1, 0x5051 << 1, 0x5054 << 1, 0x5055 << 1, 
          0x5100 << 1, 0x5101 << 1, 0x5104 << 1, 0x5105 << 1, 0x5110 << 1, 0x5111 << 1, 0x5114 << 1, 0x5115 << 1, 
          0x5140 << 1, 0x5141 << 1, 0x5144 << 1, 0x5145 << 1, 0x5150 << 1, 0x5151 << 1, 0x5154 << 1, 0x5155 << 1, 
          0x5400 << 1, 0x5401 << 1, 0x5404 << 1, 0x5405 << 1, 0x5410 << 1, 0x5411 << 1, 0x5414 << 1, 0x5415 << 1, 
          0x5440 << 1, 0x5441 << 1, 0x5444 << 1, 0x5445 << 1, 0x5450 << 1, 0x5451 << 1, 0x5454 << 1, 0x5455 << 1, 
          0x5500 << 1, 0x5501 << 1, 0x5504 << 1, 0x5505 << 1, 0x5510 << 1, 0x5511 << 1, 0x5514 << 1, 0x5515 << 1, 
          0x5540 << 1, 0x5541 << 1, 0x5544 << 1, 0x5545 << 1, 0x5550 << 1, 0x5551 << 1, 0x5554 << 1, 0x5555 << 1
        };
        static const unsigned int MortonTable256shl16[] = 
        {
          0x0000 << 16, 0x0001 << 16, 0x0004 << 16, 0x0005 << 16, 0x0010 << 16, 0x0011 << 16, 0x0014 << 16, 0x0015 << 16, 
          0x0040 << 16, 0x0041 << 16, 0x0044 << 16, 0x0045 << 16, 0x0050 << 16, 0x0051 << 16, 0x0054 << 16, 0x0055 << 16, 
          0x0100 << 16, 0x0101 << 16, 0x0104 << 16, 0x0105 << 16, 0x0110 << 16, 0x0111 << 16, 0x0114 << 16, 0x0115 << 16, 
          0x0140 << 16, 0x0141 << 16, 0x0144 << 16, 0x0145 << 16, 0x0150 << 16, 0x0151 << 16, 0x0154 << 16, 0x0155 << 16, 
          0x0400 << 16, 0x0401 << 16, 0x0404 << 16, 0x0405 << 16, 0x0410 << 16, 0x0411 << 16, 0x0414 << 16, 0x0415 << 16, 
          0x0440 << 16, 0x0441 << 16, 0x0444 << 16, 0x0445 << 16, 0x0450 << 16, 0x0451 << 16, 0x0454 << 16, 0x0455 << 16, 
          0x0500 << 16, 0x0501 << 16, 0x0504 << 16, 0x0505 << 16, 0x0510 << 16, 0x0511 << 16, 0x0514 << 16, 0x0515 << 16, 
          0x0540 << 16, 0x0541 << 16, 0x0544 << 16, 0x0545 << 16, 0x0550 << 16, 0x0551 << 16, 0x0554 << 16, 0x0555 << 16, 
          0x1000 << 16, 0x1001 << 16, 0x1004 << 16, 0x1005 << 16, 0x1010 << 16, 0x1011 << 16, 0x1014 << 16, 0x1015 << 16, 
          0x1040 << 16, 0x1041 << 16, 0x1044 << 16, 0x1045 << 16, 0x1050 << 16, 0x1051 << 16, 0x1054 << 16, 0x1055 << 16, 
          0x1100 << 16, 0x1101 << 16, 0x1104 << 16, 0x1105 << 16, 0x1110 << 16, 0x1111 << 16, 0x1114 << 16, 0x1115 << 16, 
          0x1140 << 16, 0x1141 << 16, 0x1144 << 16, 0x1145 << 16, 0x1150 << 16, 0x1151 << 16, 0x1154 << 16, 0x1155 << 16, 
          0x1400 << 16, 0x1401 << 16, 0x1404 << 16, 0x1405 << 16, 0x1410 << 16, 0x1411 << 16, 0x1414 << 16, 0x1415 << 16, 
          0x1440 << 16, 0x1441 << 16, 0x1444 << 16, 0x1445 << 16, 0x1450 << 16, 0x1451 << 16, 0x1454 << 16, 0x1455 << 16, 
          0x1500 << 16, 0x1501 << 16, 0x1504 << 16, 0x1505 << 16, 0x1510 << 16, 0x1511 << 16, 0x1514 << 16, 0x1515 << 16, 
          0x1540 << 16, 0x1541 << 16, 0x1544 << 16, 0x1545 << 16, 0x1550 << 16, 0x1551 << 16, 0x1554 << 16, 0x1555 << 16, 
          0x4000 << 16, 0x4001 << 16, 0x4004 << 16, 0x4005 << 16, 0x4010 << 16, 0x4011 << 16, 0x4014 << 16, 0x4015 << 16, 
          0x4040 << 16, 0x4041 << 16, 0x4044 << 16, 0x4045 << 16, 0x4050 << 16, 0x4051 << 16, 0x4054 << 16, 0x4055 << 16, 
          0x4100 << 16, 0x4101 << 16, 0x4104 << 16, 0x4105 << 16, 0x4110 << 16, 0x4111 << 16, 0x4114 << 16, 0x4115 << 16, 
          0x4140 << 16, 0x4141 << 16, 0x4144 << 16, 0x4145 << 16, 0x4150 << 16, 0x4151 << 16, 0x4154 << 16, 0x4155 << 16, 
          0x4400 << 16, 0x4401 << 16, 0x4404 << 16, 0x4405 << 16, 0x4410 << 16, 0x4411 << 16, 0x4414 << 16, 0x4415 << 16, 
          0x4440 << 16, 0x4441 << 16, 0x4444 << 16, 0x4445 << 16, 0x4450 << 16, 0x4451 << 16, 0x4454 << 16, 0x4455 << 16, 
          0x4500 << 16, 0x4501 << 16, 0x4504 << 16, 0x4505 << 16, 0x4510 << 16, 0x4511 << 16, 0x4514 << 16, 0x4515 << 16, 
          0x4540 << 16, 0x4541 << 16, 0x4544 << 16, 0x4545 << 16, 0x4550 << 16, 0x4551 << 16, 0x4554 << 16, 0x4555 << 16, 
          0x5000 << 16, 0x5001 << 16, 0x5004 << 16, 0x5005 << 16, 0x5010 << 16, 0x5011 << 16, 0x5014 << 16, 0x5015 << 16, 
          0x5040 << 16, 0x5041 << 16, 0x5044 << 16, 0x5045 << 16, 0x5050 << 16, 0x5051 << 16, 0x5054 << 16, 0x5055 << 16, 
          0x5100 << 16, 0x5101 << 16, 0x5104 << 16, 0x5105 << 16, 0x5110 << 16, 0x5111 << 16, 0x5114 << 16, 0x5115 << 16, 
          0x5140 << 16, 0x5141 << 16, 0x5144 << 16, 0x5145 << 16, 0x5150 << 16, 0x5151 << 16, 0x5154 << 16, 0x5155 << 16, 
          0x5400 << 16, 0x5401 << 16, 0x5404 << 16, 0x5405 << 16, 0x5410 << 16, 0x5411 << 16, 0x5414 << 16, 0x5415 << 16, 
          0x5440 << 16, 0x5441 << 16, 0x5444 << 16, 0x5445 << 16, 0x5450 << 16, 0x5451 << 16, 0x5454 << 16, 0x5455 << 16, 
          0x5500 << 16, 0x5501 << 16, 0x5504 << 16, 0x5505 << 16, 0x5510 << 16, 0x5511 << 16, 0x5514 << 16, 0x5515 << 16, 
          0x5540 << 16, 0x5541 << 16, 0x5544 << 16, 0x5545 << 16, 0x5550 << 16, 0x5551 << 16, 0x5554 << 16, 0x5555 << 16
        };
        static const unsigned int MortonTable256shl17[] = 
        {
          0x0000 << 17, 0x0001 << 17, 0x0004 << 17, 0x0005 << 17, 0x0010 << 17, 0x0011 << 17, 0x0014 << 17, 0x0015 << 17, 
          0x0040 << 17, 0x0041 << 17, 0x0044 << 17, 0x0045 << 17, 0x0050 << 17, 0x0051 << 17, 0x0054 << 17, 0x0055 << 17, 
          0x0100 << 17, 0x0101 << 17, 0x0104 << 17, 0x0105 << 17, 0x0110 << 17, 0x0111 << 17, 0x0114 << 17, 0x0115 << 17, 
          0x0140 << 17, 0x0141 << 17, 0x0144 << 17, 0x0145 << 17, 0x0150 << 17, 0x0151 << 17, 0x0154 << 17, 0x0155 << 17, 
          0x0400 << 17, 0x0401 << 17, 0x0404 << 17, 0x0405 << 17, 0x0410 << 17, 0x0411 << 17, 0x0414 << 17, 0x0415 << 17, 
          0x0440 << 17, 0x0441 << 17, 0x0444 << 17, 0x0445 << 17, 0x0450 << 17, 0x0451 << 17, 0x0454 << 17, 0x0455 << 17, 
          0x0500 << 17, 0x0501 << 17, 0x0504 << 17, 0x0505 << 17, 0x0510 << 17, 0x0511 << 17, 0x0514 << 17, 0x0515 << 17, 
          0x0540 << 17, 0x0541 << 17, 0x0544 << 17, 0x0545 << 17, 0x0550 << 17, 0x0551 << 17, 0x0554 << 17, 0x0555 << 17, 
          0x1000 << 17, 0x1001 << 17, 0x1004 << 17, 0x1005 << 17, 0x1010 << 17, 0x1011 << 17, 0x1014 << 17, 0x1015 << 17, 
          0x1040 << 17, 0x1041 << 17, 0x1044 << 17, 0x1045 << 17, 0x1050 << 17, 0x1051 << 17, 0x1054 << 17, 0x1055 << 17, 
          0x1100 << 17, 0x1101 << 17, 0x1104 << 17, 0x1105 << 17, 0x1110 << 17, 0x1111 << 17, 0x1114 << 17, 0x1115 << 17, 
          0x1140 << 17, 0x1141 << 17, 0x1144 << 17, 0x1145 << 17, 0x1150 << 17, 0x1151 << 17, 0x1154 << 17, 0x1155 << 17, 
          0x1400 << 17, 0x1401 << 17, 0x1404 << 17, 0x1405 << 17, 0x1410 << 17, 0x1411 << 17, 0x1414 << 17, 0x1415 << 17, 
          0x1440 << 17, 0x1441 << 17, 0x1444 << 17, 0x1445 << 17, 0x1450 << 17, 0x1451 << 17, 0x1454 << 17, 0x1455 << 17, 
          0x1500 << 17, 0x1501 << 17, 0x1504 << 17, 0x1505 << 17, 0x1510 << 17, 0x1511 << 17, 0x1514 << 17, 0x1515 << 17, 
          0x1540 << 17, 0x1541 << 17, 0x1544 << 17, 0x1545 << 17, 0x1550 << 17, 0x1551 << 17, 0x1554 << 17, 0x1555 << 17, 
          0x4000 << 17, 0x4001 << 17, 0x4004 << 17, 0x4005 << 17, 0x4010 << 17, 0x4011 << 17, 0x4014 << 17, 0x4015 << 17, 
          0x4040 << 17, 0x4041 << 17, 0x4044 << 17, 0x4045 << 17, 0x4050 << 17, 0x4051 << 17, 0x4054 << 17, 0x4055 << 17, 
          0x4100 << 17, 0x4101 << 17, 0x4104 << 17, 0x4105 << 17, 0x4110 << 17, 0x4111 << 17, 0x4114 << 17, 0x4115 << 17, 
          0x4140 << 17, 0x4141 << 17, 0x4144 << 17, 0x4145 << 17, 0x4150 << 17, 0x4151 << 17, 0x4154 << 17, 0x4155 << 17, 
          0x4400 << 17, 0x4401 << 17, 0x4404 << 17, 0x4405 << 17, 0x4410 << 17, 0x4411 << 17, 0x4414 << 17, 0x4415 << 17, 
          0x4440 << 17, 0x4441 << 17, 0x4444 << 17, 0x4445 << 17, 0x4450 << 17, 0x4451 << 17, 0x4454 << 17, 0x4455 << 17, 
          0x4500 << 17, 0x4501 << 17, 0x4504 << 17, 0x4505 << 17, 0x4510 << 17, 0x4511 << 17, 0x4514 << 17, 0x4515 << 17, 
          0x4540 << 17, 0x4541 << 17, 0x4544 << 17, 0x4545 << 17, 0x4550 << 17, 0x4551 << 17, 0x4554 << 17, 0x4555 << 17, 
          0x5000 << 17, 0x5001 << 17, 0x5004 << 17, 0x5005 << 17, 0x5010 << 17, 0x5011 << 17, 0x5014 << 17, 0x5015 << 17, 
          0x5040 << 17, 0x5041 << 17, 0x5044 << 17, 0x5045 << 17, 0x5050 << 17, 0x5051 << 17, 0x5054 << 17, 0x5055 << 17, 
          0x5100 << 17, 0x5101 << 17, 0x5104 << 17, 0x5105 << 17, 0x5110 << 17, 0x5111 << 17, 0x5114 << 17, 0x5115 << 17, 
          0x5140 << 17, 0x5141 << 17, 0x5144 << 17, 0x5145 << 17, 0x5150 << 17, 0x5151 << 17, 0x5154 << 17, 0x5155 << 17, 
          0x5400 << 17, 0x5401 << 17, 0x5404 << 17, 0x5405 << 17, 0x5410 << 17, 0x5411 << 17, 0x5414 << 17, 0x5415 << 17, 
          0x5440 << 17, 0x5441 << 17, 0x5444 << 17, 0x5445 << 17, 0x5450 << 17, 0x5451 << 17, 0x5454 << 17, 0x5455 << 17, 
          0x5500 << 17, 0x5501 << 17, 0x5504 << 17, 0x5505 << 17, 0x5510 << 17, 0x5511 << 17, 0x5514 << 17, 0x5515 << 17, 
          0x5540 << 17, 0x5541 << 17, 0x5544 << 17, 0x5545 << 17, 0x5550 << 17, 0x5551 << 17, 0x5554 << 17, 0x5555 << 17
        };
}

template <int size> struct MakeMortonIndex {
private:

public:
        size_t operator () (unsigned short x, unsigned short y) const __attribute__((hot,const)){
                //static const unsigned int *MortonTable256_ofs = MortonTable256 - 
                using Morton::MortonTable256;
                return MortonTable256[y >> 8] << 17 | 
                       MortonTable256[x >> 8] << 16 |
                       MortonTable256[y & 0xFF] << 1 | 
                       MortonTable256[x & 0xFF];

        }
        
        size_t computeX (unsigned short x) const __attribute__((hot,const)){
                using Morton::MortonTable256;
                return MortonTable256[x >> 8] << 16 |                       
                       MortonTable256[x & 0xFF];
        }
        
        size_t computeY (unsigned short y) const __attribute__((hot,const)){
                using Morton::MortonTable256;
                return MortonTable256[y >> 8] << 17 |                        
                       MortonTable256[y & 0xFF] << 1;
        }
        
        size_t mergeXY (size_t X, size_t Y) const __attribute__((hot,const)) {
                return X | Y;
        }
};


inline float f2_to_hemisphere (float fx, float fy) {
        const float fz_sq = 1.0 - fx*fx - fy*fy;
        if (fz_sq < 0.0) {
                return -1;
        }
        return sqrt (fz_sq);
}

//-----
#define StdMakeIndex MakeMortonIndex


#define SSE
#include <sse.hh>


#include "preetham.hh"



///////////////////////////////////////////////////////////////////////////////
// Renderer
///////////////////////////////////////////////////////////////////////////////
class Renderer {
public:
        Renderer (unsigned int width_, unsigned int height_);
        ~Renderer ();
        void render (float,float,float);
        bool isInitialized () const;

        
        struct ray_t {
                float position [4];
                float direction [4];                
                
                ray_t (const float *position_, 
                        const float *direction_
                ) {
                        for (int i=0; i<4; ++i) {
                                position [i] = position_ [i];
                                direction [i] = direction_ [i];
                        }
                }
                
                ray_t() {}
        private:
        };

private:   

        void renderTask () const;


        // ...
        float currentPosition [4];
        float currentYaw;

        const unsigned int width,height;
        SDL_Surface* screen;

        
        //ray_t generateRay (int x, int y) const ;
        //float intersect (ray_t rayf, float depth = 0.0f) const ;

        const float* skylight (float x, float y, float z) const ;

        struct Vertex1 {
                float height;
                float normal [3];
        };
        Vertex1 fun (float x, float y) const ;        

        #ifdef SSE
        typedef grind::sse::Float Float;
        typedef grind::sse::Vector3d Vector3d;
        typedef grind::sse::Ray Ray;
        typedef grind::sse::Mask Mask;

        struct Vertex {                
                Float height;
                Vector3d normal;
                
                Vertex () {}
                Vertex (const Float &height_, const Vector3d &normal_) 
                : height (height_), normal (normal_) {}
        };

        struct Intersection {
                Float depth;
                Vertex vertex;
                
                Intersection () {}
                Intersection (const Float &depth_, const Vertex &vertex_) 
                : depth (depth_), vertex (vertex_) {}
        };
        Intersection intersect (Ray ray, Float depth, Mask active) const __attribute__((hot));//,noinline));
        Vertex fun (Float x, Float y) const __attribute__((hot));
        Ray generateRay4 (int x, int y) const __attribute__((hot));
        void generateRayDirection (float dir[3], int x, int y) const __attribute__((hot));
        #endif
        
        
        // heightmap
        enum { heightfield_size = 1024 };
        Vertex1 *heightfield;
        void initHeightfield () ;
        StdMakeIndex<heightfield_size> makeHeightfieldIndex;
        const float heightfield_scale;
        
        // skymap
        enum { skymap_size = 1024 };
        StdMakeIndex<skymap_size> makeSkymapIndex;
        float skymap_fsize;
        struct rgbx { float r, g, b, x; };
        float (*skymap) [4];
        Vector3d sunDirection;
        
        void initSkymap();
};

#endif // RENDERER_HH_INCLUDED_20090503