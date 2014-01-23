-- (C) 1992-2013 Altera Corporation. All rights reserved.                         
-- Your use of Altera Corporation's design tools, logic functions and other       
-- software and tools, and its AMPP partner logic functions, and any output       
-- files any of the foregoing (including device programming or simulation         
-- files), and any associated documentation or information are expressly subject  
-- to the terms and conditions of the Altera Program License Subscription         
-- Agreement, Altera MegaCore Function License Agreement, or other applicable     
-- license agreement, including, without limitation, that your use is for the     
-- sole purpose of programming logic devices manufactured by Altera and sold by   
-- Altera or its authorized distributors.  Please refer to the applicable         
-- agreement for further details.                                                 
    


LIBRARY ieee;
USE ieee.std_logic_1164.all;
USE ieee.std_logic_unsigned.all;
USE ieee.std_logic_arith.all;

--***************************************************
--***                                             ***
--***   FLOATING POINT CORE LIBRARY               ***
--***                                             ***
--***   FP_INVSQR_LUT1.VHD                        ***
--***                                             ***
--***   Function: Look Up Table - Inverse Root    ***
--***                                             ***
--***   Generated by MATLAB Utility               ***
--***                                             ***
--***   31/01/08 ML                               ***
--***                                             ***
--***   (c) 2008 Altera Corporation               ***
--***                                             ***
--***   Change History                            ***
--***                                             ***
--***                                             ***
--***                                             ***
--***                                             ***
--***                                             ***
--***************************************************

ENTITY fp_invsqr_lut1 IS
PORT (
      add : IN STD_LOGIC_VECTOR (9 DOWNTO 1);
		data : OUT STD_LOGIC_VECTOR (11 DOWNTO 1)
);
END fp_invsqr_lut1;

ARCHITECTURE rtl OF fp_invsqr_lut1 IS

BEGIN

  pca: PROCESS (add)
  BEGIN
    CASE add IS
      WHEN "000000000" => data <= conv_std_logic_vector(1023,11);
      WHEN "000000001" => data <= conv_std_logic_vector(1020,11);
      WHEN "000000010" => data <= conv_std_logic_vector(1017,11);
      WHEN "000000011" => data <= conv_std_logic_vector(1014,11);
      WHEN "000000100" => data <= conv_std_logic_vector(1011,11);
      WHEN "000000101" => data <= conv_std_logic_vector(1008,11);
      WHEN "000000110" => data <= conv_std_logic_vector(1005,11);
      WHEN "000000111" => data <= conv_std_logic_vector(1002,11);
      WHEN "000001000" => data <= conv_std_logic_vector(999,11);
      WHEN "000001001" => data <= conv_std_logic_vector(996,11);
      WHEN "000001010" => data <= conv_std_logic_vector(993,11);
      WHEN "000001011" => data <= conv_std_logic_vector(990,11);
      WHEN "000001100" => data <= conv_std_logic_vector(988,11);
      WHEN "000001101" => data <= conv_std_logic_vector(985,11);
      WHEN "000001110" => data <= conv_std_logic_vector(982,11);
      WHEN "000001111" => data <= conv_std_logic_vector(979,11);
      WHEN "000010000" => data <= conv_std_logic_vector(976,11);
      WHEN "000010001" => data <= conv_std_logic_vector(974,11);
      WHEN "000010010" => data <= conv_std_logic_vector(971,11);
      WHEN "000010011" => data <= conv_std_logic_vector(968,11);
      WHEN "000010100" => data <= conv_std_logic_vector(965,11);
      WHEN "000010101" => data <= conv_std_logic_vector(963,11);
      WHEN "000010110" => data <= conv_std_logic_vector(960,11);
      WHEN "000010111" => data <= conv_std_logic_vector(957,11);
      WHEN "000011000" => data <= conv_std_logic_vector(955,11);
      WHEN "000011001" => data <= conv_std_logic_vector(952,11);
      WHEN "000011010" => data <= conv_std_logic_vector(949,11);
      WHEN "000011011" => data <= conv_std_logic_vector(947,11);
      WHEN "000011100" => data <= conv_std_logic_vector(944,11);
      WHEN "000011101" => data <= conv_std_logic_vector(941,11);
      WHEN "000011110" => data <= conv_std_logic_vector(939,11);
      WHEN "000011111" => data <= conv_std_logic_vector(936,11);
      WHEN "000100000" => data <= conv_std_logic_vector(934,11);
      WHEN "000100001" => data <= conv_std_logic_vector(931,11);
      WHEN "000100010" => data <= conv_std_logic_vector(929,11);
      WHEN "000100011" => data <= conv_std_logic_vector(926,11);
      WHEN "000100100" => data <= conv_std_logic_vector(924,11);
      WHEN "000100101" => data <= conv_std_logic_vector(921,11);
      WHEN "000100110" => data <= conv_std_logic_vector(918,11);
      WHEN "000100111" => data <= conv_std_logic_vector(916,11);
      WHEN "000101000" => data <= conv_std_logic_vector(913,11);
      WHEN "000101001" => data <= conv_std_logic_vector(911,11);
      WHEN "000101010" => data <= conv_std_logic_vector(909,11);
      WHEN "000101011" => data <= conv_std_logic_vector(906,11);
      WHEN "000101100" => data <= conv_std_logic_vector(904,11);
      WHEN "000101101" => data <= conv_std_logic_vector(901,11);
      WHEN "000101110" => data <= conv_std_logic_vector(899,11);
      WHEN "000101111" => data <= conv_std_logic_vector(896,11);
      WHEN "000110000" => data <= conv_std_logic_vector(894,11);
      WHEN "000110001" => data <= conv_std_logic_vector(892,11);
      WHEN "000110010" => data <= conv_std_logic_vector(889,11);
      WHEN "000110011" => data <= conv_std_logic_vector(887,11);
      WHEN "000110100" => data <= conv_std_logic_vector(885,11);
      WHEN "000110101" => data <= conv_std_logic_vector(882,11);
      WHEN "000110110" => data <= conv_std_logic_vector(880,11);
      WHEN "000110111" => data <= conv_std_logic_vector(878,11);
      WHEN "000111000" => data <= conv_std_logic_vector(875,11);
      WHEN "000111001" => data <= conv_std_logic_vector(873,11);
      WHEN "000111010" => data <= conv_std_logic_vector(871,11);
      WHEN "000111011" => data <= conv_std_logic_vector(868,11);
      WHEN "000111100" => data <= conv_std_logic_vector(866,11);
      WHEN "000111101" => data <= conv_std_logic_vector(864,11);
      WHEN "000111110" => data <= conv_std_logic_vector(862,11);
      WHEN "000111111" => data <= conv_std_logic_vector(859,11);
      WHEN "001000000" => data <= conv_std_logic_vector(857,11);
      WHEN "001000001" => data <= conv_std_logic_vector(855,11);
      WHEN "001000010" => data <= conv_std_logic_vector(853,11);
      WHEN "001000011" => data <= conv_std_logic_vector(850,11);
      WHEN "001000100" => data <= conv_std_logic_vector(848,11);
      WHEN "001000101" => data <= conv_std_logic_vector(846,11);
      WHEN "001000110" => data <= conv_std_logic_vector(844,11);
      WHEN "001000111" => data <= conv_std_logic_vector(842,11);
      WHEN "001001000" => data <= conv_std_logic_vector(840,11);
      WHEN "001001001" => data <= conv_std_logic_vector(837,11);
      WHEN "001001010" => data <= conv_std_logic_vector(835,11);
      WHEN "001001011" => data <= conv_std_logic_vector(833,11);
      WHEN "001001100" => data <= conv_std_logic_vector(831,11);
      WHEN "001001101" => data <= conv_std_logic_vector(829,11);
      WHEN "001001110" => data <= conv_std_logic_vector(827,11);
      WHEN "001001111" => data <= conv_std_logic_vector(825,11);
      WHEN "001010000" => data <= conv_std_logic_vector(823,11);
      WHEN "001010001" => data <= conv_std_logic_vector(820,11);
      WHEN "001010010" => data <= conv_std_logic_vector(818,11);
      WHEN "001010011" => data <= conv_std_logic_vector(816,11);
      WHEN "001010100" => data <= conv_std_logic_vector(814,11);
      WHEN "001010101" => data <= conv_std_logic_vector(812,11);
      WHEN "001010110" => data <= conv_std_logic_vector(810,11);
      WHEN "001010111" => data <= conv_std_logic_vector(808,11);
      WHEN "001011000" => data <= conv_std_logic_vector(806,11);
      WHEN "001011001" => data <= conv_std_logic_vector(804,11);
      WHEN "001011010" => data <= conv_std_logic_vector(802,11);
      WHEN "001011011" => data <= conv_std_logic_vector(800,11);
      WHEN "001011100" => data <= conv_std_logic_vector(798,11);
      WHEN "001011101" => data <= conv_std_logic_vector(796,11);
      WHEN "001011110" => data <= conv_std_logic_vector(794,11);
      WHEN "001011111" => data <= conv_std_logic_vector(792,11);
      WHEN "001100000" => data <= conv_std_logic_vector(790,11);
      WHEN "001100001" => data <= conv_std_logic_vector(788,11);
      WHEN "001100010" => data <= conv_std_logic_vector(786,11);
      WHEN "001100011" => data <= conv_std_logic_vector(785,11);
      WHEN "001100100" => data <= conv_std_logic_vector(783,11);
      WHEN "001100101" => data <= conv_std_logic_vector(781,11);
      WHEN "001100110" => data <= conv_std_logic_vector(779,11);
      WHEN "001100111" => data <= conv_std_logic_vector(777,11);
      WHEN "001101000" => data <= conv_std_logic_vector(775,11);
      WHEN "001101001" => data <= conv_std_logic_vector(773,11);
      WHEN "001101010" => data <= conv_std_logic_vector(771,11);
      WHEN "001101011" => data <= conv_std_logic_vector(769,11);
      WHEN "001101100" => data <= conv_std_logic_vector(768,11);
      WHEN "001101101" => data <= conv_std_logic_vector(766,11);
      WHEN "001101110" => data <= conv_std_logic_vector(764,11);
      WHEN "001101111" => data <= conv_std_logic_vector(762,11);
      WHEN "001110000" => data <= conv_std_logic_vector(760,11);
      WHEN "001110001" => data <= conv_std_logic_vector(758,11);
      WHEN "001110010" => data <= conv_std_logic_vector(757,11);
      WHEN "001110011" => data <= conv_std_logic_vector(755,11);
      WHEN "001110100" => data <= conv_std_logic_vector(753,11);
      WHEN "001110101" => data <= conv_std_logic_vector(751,11);
      WHEN "001110110" => data <= conv_std_logic_vector(749,11);
      WHEN "001110111" => data <= conv_std_logic_vector(748,11);
      WHEN "001111000" => data <= conv_std_logic_vector(746,11);
      WHEN "001111001" => data <= conv_std_logic_vector(744,11);
      WHEN "001111010" => data <= conv_std_logic_vector(742,11);
      WHEN "001111011" => data <= conv_std_logic_vector(741,11);
      WHEN "001111100" => data <= conv_std_logic_vector(739,11);
      WHEN "001111101" => data <= conv_std_logic_vector(737,11);
      WHEN "001111110" => data <= conv_std_logic_vector(735,11);
      WHEN "001111111" => data <= conv_std_logic_vector(734,11);
      WHEN "010000000" => data <= conv_std_logic_vector(732,11);
      WHEN "010000001" => data <= conv_std_logic_vector(730,11);
      WHEN "010000010" => data <= conv_std_logic_vector(728,11);
      WHEN "010000011" => data <= conv_std_logic_vector(727,11);
      WHEN "010000100" => data <= conv_std_logic_vector(725,11);
      WHEN "010000101" => data <= conv_std_logic_vector(723,11);
      WHEN "010000110" => data <= conv_std_logic_vector(722,11);
      WHEN "010000111" => data <= conv_std_logic_vector(720,11);
      WHEN "010001000" => data <= conv_std_logic_vector(718,11);
      WHEN "010001001" => data <= conv_std_logic_vector(717,11);
      WHEN "010001010" => data <= conv_std_logic_vector(715,11);
      WHEN "010001011" => data <= conv_std_logic_vector(713,11);
      WHEN "010001100" => data <= conv_std_logic_vector(712,11);
      WHEN "010001101" => data <= conv_std_logic_vector(710,11);
      WHEN "010001110" => data <= conv_std_logic_vector(709,11);
      WHEN "010001111" => data <= conv_std_logic_vector(707,11);
      WHEN "010010000" => data <= conv_std_logic_vector(705,11);
      WHEN "010010001" => data <= conv_std_logic_vector(704,11);
      WHEN "010010010" => data <= conv_std_logic_vector(702,11);
      WHEN "010010011" => data <= conv_std_logic_vector(700,11);
      WHEN "010010100" => data <= conv_std_logic_vector(699,11);
      WHEN "010010101" => data <= conv_std_logic_vector(697,11);
      WHEN "010010110" => data <= conv_std_logic_vector(696,11);
      WHEN "010010111" => data <= conv_std_logic_vector(694,11);
      WHEN "010011000" => data <= conv_std_logic_vector(693,11);
      WHEN "010011001" => data <= conv_std_logic_vector(691,11);
      WHEN "010011010" => data <= conv_std_logic_vector(689,11);
      WHEN "010011011" => data <= conv_std_logic_vector(688,11);
      WHEN "010011100" => data <= conv_std_logic_vector(686,11);
      WHEN "010011101" => data <= conv_std_logic_vector(685,11);
      WHEN "010011110" => data <= conv_std_logic_vector(683,11);
      WHEN "010011111" => data <= conv_std_logic_vector(682,11);
      WHEN "010100000" => data <= conv_std_logic_vector(680,11);
      WHEN "010100001" => data <= conv_std_logic_vector(679,11);
      WHEN "010100010" => data <= conv_std_logic_vector(677,11);
      WHEN "010100011" => data <= conv_std_logic_vector(676,11);
      WHEN "010100100" => data <= conv_std_logic_vector(674,11);
      WHEN "010100101" => data <= conv_std_logic_vector(673,11);
      WHEN "010100110" => data <= conv_std_logic_vector(671,11);
      WHEN "010100111" => data <= conv_std_logic_vector(670,11);
      WHEN "010101000" => data <= conv_std_logic_vector(668,11);
      WHEN "010101001" => data <= conv_std_logic_vector(667,11);
      WHEN "010101010" => data <= conv_std_logic_vector(665,11);
      WHEN "010101011" => data <= conv_std_logic_vector(664,11);
      WHEN "010101100" => data <= conv_std_logic_vector(662,11);
      WHEN "010101101" => data <= conv_std_logic_vector(661,11);
      WHEN "010101110" => data <= conv_std_logic_vector(660,11);
      WHEN "010101111" => data <= conv_std_logic_vector(658,11);
      WHEN "010110000" => data <= conv_std_logic_vector(657,11);
      WHEN "010110001" => data <= conv_std_logic_vector(655,11);
      WHEN "010110010" => data <= conv_std_logic_vector(654,11);
      WHEN "010110011" => data <= conv_std_logic_vector(652,11);
      WHEN "010110100" => data <= conv_std_logic_vector(651,11);
      WHEN "010110101" => data <= conv_std_logic_vector(650,11);
      WHEN "010110110" => data <= conv_std_logic_vector(648,11);
      WHEN "010110111" => data <= conv_std_logic_vector(647,11);
      WHEN "010111000" => data <= conv_std_logic_vector(645,11);
      WHEN "010111001" => data <= conv_std_logic_vector(644,11);
      WHEN "010111010" => data <= conv_std_logic_vector(643,11);
      WHEN "010111011" => data <= conv_std_logic_vector(641,11);
      WHEN "010111100" => data <= conv_std_logic_vector(640,11);
      WHEN "010111101" => data <= conv_std_logic_vector(639,11);
      WHEN "010111110" => data <= conv_std_logic_vector(637,11);
      WHEN "010111111" => data <= conv_std_logic_vector(636,11);
      WHEN "011000000" => data <= conv_std_logic_vector(634,11);
      WHEN "011000001" => data <= conv_std_logic_vector(633,11);
      WHEN "011000010" => data <= conv_std_logic_vector(632,11);
      WHEN "011000011" => data <= conv_std_logic_vector(630,11);
      WHEN "011000100" => data <= conv_std_logic_vector(629,11);
      WHEN "011000101" => data <= conv_std_logic_vector(628,11);
      WHEN "011000110" => data <= conv_std_logic_vector(626,11);
      WHEN "011000111" => data <= conv_std_logic_vector(625,11);
      WHEN "011001000" => data <= conv_std_logic_vector(624,11);
      WHEN "011001001" => data <= conv_std_logic_vector(622,11);
      WHEN "011001010" => data <= conv_std_logic_vector(621,11);
      WHEN "011001011" => data <= conv_std_logic_vector(620,11);
      WHEN "011001100" => data <= conv_std_logic_vector(619,11);
      WHEN "011001101" => data <= conv_std_logic_vector(617,11);
      WHEN "011001110" => data <= conv_std_logic_vector(616,11);
      WHEN "011001111" => data <= conv_std_logic_vector(615,11);
      WHEN "011010000" => data <= conv_std_logic_vector(613,11);
      WHEN "011010001" => data <= conv_std_logic_vector(612,11);
      WHEN "011010010" => data <= conv_std_logic_vector(611,11);
      WHEN "011010011" => data <= conv_std_logic_vector(610,11);
      WHEN "011010100" => data <= conv_std_logic_vector(608,11);
      WHEN "011010101" => data <= conv_std_logic_vector(607,11);
      WHEN "011010110" => data <= conv_std_logic_vector(606,11);
      WHEN "011010111" => data <= conv_std_logic_vector(605,11);
      WHEN "011011000" => data <= conv_std_logic_vector(603,11);
      WHEN "011011001" => data <= conv_std_logic_vector(602,11);
      WHEN "011011010" => data <= conv_std_logic_vector(601,11);
      WHEN "011011011" => data <= conv_std_logic_vector(600,11);
      WHEN "011011100" => data <= conv_std_logic_vector(598,11);
      WHEN "011011101" => data <= conv_std_logic_vector(597,11);
      WHEN "011011110" => data <= conv_std_logic_vector(596,11);
      WHEN "011011111" => data <= conv_std_logic_vector(595,11);
      WHEN "011100000" => data <= conv_std_logic_vector(594,11);
      WHEN "011100001" => data <= conv_std_logic_vector(592,11);
      WHEN "011100010" => data <= conv_std_logic_vector(591,11);
      WHEN "011100011" => data <= conv_std_logic_vector(590,11);
      WHEN "011100100" => data <= conv_std_logic_vector(589,11);
      WHEN "011100101" => data <= conv_std_logic_vector(588,11);
      WHEN "011100110" => data <= conv_std_logic_vector(586,11);
      WHEN "011100111" => data <= conv_std_logic_vector(585,11);
      WHEN "011101000" => data <= conv_std_logic_vector(584,11);
      WHEN "011101001" => data <= conv_std_logic_vector(583,11);
      WHEN "011101010" => data <= conv_std_logic_vector(582,11);
      WHEN "011101011" => data <= conv_std_logic_vector(580,11);
      WHEN "011101100" => data <= conv_std_logic_vector(579,11);
      WHEN "011101101" => data <= conv_std_logic_vector(578,11);
      WHEN "011101110" => data <= conv_std_logic_vector(577,11);
      WHEN "011101111" => data <= conv_std_logic_vector(576,11);
      WHEN "011110000" => data <= conv_std_logic_vector(575,11);
      WHEN "011110001" => data <= conv_std_logic_vector(574,11);
      WHEN "011110010" => data <= conv_std_logic_vector(572,11);
      WHEN "011110011" => data <= conv_std_logic_vector(571,11);
      WHEN "011110100" => data <= conv_std_logic_vector(570,11);
      WHEN "011110101" => data <= conv_std_logic_vector(569,11);
      WHEN "011110110" => data <= conv_std_logic_vector(568,11);
      WHEN "011110111" => data <= conv_std_logic_vector(567,11);
      WHEN "011111000" => data <= conv_std_logic_vector(566,11);
      WHEN "011111001" => data <= conv_std_logic_vector(565,11);
      WHEN "011111010" => data <= conv_std_logic_vector(563,11);
      WHEN "011111011" => data <= conv_std_logic_vector(562,11);
      WHEN "011111100" => data <= conv_std_logic_vector(561,11);
      WHEN "011111101" => data <= conv_std_logic_vector(560,11);
      WHEN "011111110" => data <= conv_std_logic_vector(559,11);
      WHEN "011111111" => data <= conv_std_logic_vector(558,11);
      WHEN "100000000" => data <= conv_std_logic_vector(557,11);
      WHEN "100000001" => data <= conv_std_logic_vector(556,11);
      WHEN "100000010" => data <= conv_std_logic_vector(555,11);
      WHEN "100000011" => data <= conv_std_logic_vector(554,11);
      WHEN "100000100" => data <= conv_std_logic_vector(553,11);
      WHEN "100000101" => data <= conv_std_logic_vector(551,11);
      WHEN "100000110" => data <= conv_std_logic_vector(550,11);
      WHEN "100000111" => data <= conv_std_logic_vector(549,11);
      WHEN "100001000" => data <= conv_std_logic_vector(548,11);
      WHEN "100001001" => data <= conv_std_logic_vector(547,11);
      WHEN "100001010" => data <= conv_std_logic_vector(546,11);
      WHEN "100001011" => data <= conv_std_logic_vector(545,11);
      WHEN "100001100" => data <= conv_std_logic_vector(544,11);
      WHEN "100001101" => data <= conv_std_logic_vector(543,11);
      WHEN "100001110" => data <= conv_std_logic_vector(542,11);
      WHEN "100001111" => data <= conv_std_logic_vector(541,11);
      WHEN "100010000" => data <= conv_std_logic_vector(540,11);
      WHEN "100010001" => data <= conv_std_logic_vector(539,11);
      WHEN "100010010" => data <= conv_std_logic_vector(538,11);
      WHEN "100010011" => data <= conv_std_logic_vector(537,11);
      WHEN "100010100" => data <= conv_std_logic_vector(536,11);
      WHEN "100010101" => data <= conv_std_logic_vector(535,11);
      WHEN "100010110" => data <= conv_std_logic_vector(534,11);
      WHEN "100010111" => data <= conv_std_logic_vector(533,11);
      WHEN "100011000" => data <= conv_std_logic_vector(532,11);
      WHEN "100011001" => data <= conv_std_logic_vector(531,11);
      WHEN "100011010" => data <= conv_std_logic_vector(530,11);
      WHEN "100011011" => data <= conv_std_logic_vector(529,11);
      WHEN "100011100" => data <= conv_std_logic_vector(528,11);
      WHEN "100011101" => data <= conv_std_logic_vector(527,11);
      WHEN "100011110" => data <= conv_std_logic_vector(526,11);
      WHEN "100011111" => data <= conv_std_logic_vector(525,11);
      WHEN "100100000" => data <= conv_std_logic_vector(524,11);
      WHEN "100100001" => data <= conv_std_logic_vector(523,11);
      WHEN "100100010" => data <= conv_std_logic_vector(522,11);
      WHEN "100100011" => data <= conv_std_logic_vector(521,11);
      WHEN "100100100" => data <= conv_std_logic_vector(520,11);
      WHEN "100100101" => data <= conv_std_logic_vector(519,11);
      WHEN "100100110" => data <= conv_std_logic_vector(518,11);
      WHEN "100100111" => data <= conv_std_logic_vector(517,11);
      WHEN "100101000" => data <= conv_std_logic_vector(516,11);
      WHEN "100101001" => data <= conv_std_logic_vector(515,11);
      WHEN "100101010" => data <= conv_std_logic_vector(514,11);
      WHEN "100101011" => data <= conv_std_logic_vector(513,11);
      WHEN "100101100" => data <= conv_std_logic_vector(512,11);
      WHEN "100101101" => data <= conv_std_logic_vector(511,11);
      WHEN "100101110" => data <= conv_std_logic_vector(510,11);
      WHEN "100101111" => data <= conv_std_logic_vector(509,11);
      WHEN "100110000" => data <= conv_std_logic_vector(508,11);
      WHEN "100110001" => data <= conv_std_logic_vector(508,11);
      WHEN "100110010" => data <= conv_std_logic_vector(507,11);
      WHEN "100110011" => data <= conv_std_logic_vector(506,11);
      WHEN "100110100" => data <= conv_std_logic_vector(505,11);
      WHEN "100110101" => data <= conv_std_logic_vector(504,11);
      WHEN "100110110" => data <= conv_std_logic_vector(503,11);
      WHEN "100110111" => data <= conv_std_logic_vector(502,11);
      WHEN "100111000" => data <= conv_std_logic_vector(501,11);
      WHEN "100111001" => data <= conv_std_logic_vector(500,11);
      WHEN "100111010" => data <= conv_std_logic_vector(499,11);
      WHEN "100111011" => data <= conv_std_logic_vector(498,11);
      WHEN "100111100" => data <= conv_std_logic_vector(497,11);
      WHEN "100111101" => data <= conv_std_logic_vector(497,11);
      WHEN "100111110" => data <= conv_std_logic_vector(496,11);
      WHEN "100111111" => data <= conv_std_logic_vector(495,11);
      WHEN "101000000" => data <= conv_std_logic_vector(494,11);
      WHEN "101000001" => data <= conv_std_logic_vector(493,11);
      WHEN "101000010" => data <= conv_std_logic_vector(492,11);
      WHEN "101000011" => data <= conv_std_logic_vector(491,11);
      WHEN "101000100" => data <= conv_std_logic_vector(490,11);
      WHEN "101000101" => data <= conv_std_logic_vector(489,11);
      WHEN "101000110" => data <= conv_std_logic_vector(489,11);
      WHEN "101000111" => data <= conv_std_logic_vector(488,11);
      WHEN "101001000" => data <= conv_std_logic_vector(487,11);
      WHEN "101001001" => data <= conv_std_logic_vector(486,11);
      WHEN "101001010" => data <= conv_std_logic_vector(485,11);
      WHEN "101001011" => data <= conv_std_logic_vector(484,11);
      WHEN "101001100" => data <= conv_std_logic_vector(483,11);
      WHEN "101001101" => data <= conv_std_logic_vector(483,11);
      WHEN "101001110" => data <= conv_std_logic_vector(482,11);
      WHEN "101001111" => data <= conv_std_logic_vector(481,11);
      WHEN "101010000" => data <= conv_std_logic_vector(480,11);
      WHEN "101010001" => data <= conv_std_logic_vector(479,11);
      WHEN "101010010" => data <= conv_std_logic_vector(478,11);
      WHEN "101010011" => data <= conv_std_logic_vector(477,11);
      WHEN "101010100" => data <= conv_std_logic_vector(477,11);
      WHEN "101010101" => data <= conv_std_logic_vector(476,11);
      WHEN "101010110" => data <= conv_std_logic_vector(475,11);
      WHEN "101010111" => data <= conv_std_logic_vector(474,11);
      WHEN "101011000" => data <= conv_std_logic_vector(473,11);
      WHEN "101011001" => data <= conv_std_logic_vector(472,11);
      WHEN "101011010" => data <= conv_std_logic_vector(472,11);
      WHEN "101011011" => data <= conv_std_logic_vector(471,11);
      WHEN "101011100" => data <= conv_std_logic_vector(470,11);
      WHEN "101011101" => data <= conv_std_logic_vector(469,11);
      WHEN "101011110" => data <= conv_std_logic_vector(468,11);
      WHEN "101011111" => data <= conv_std_logic_vector(468,11);
      WHEN "101100000" => data <= conv_std_logic_vector(467,11);
      WHEN "101100001" => data <= conv_std_logic_vector(466,11);
      WHEN "101100010" => data <= conv_std_logic_vector(465,11);
      WHEN "101100011" => data <= conv_std_logic_vector(464,11);
      WHEN "101100100" => data <= conv_std_logic_vector(464,11);
      WHEN "101100101" => data <= conv_std_logic_vector(463,11);
      WHEN "101100110" => data <= conv_std_logic_vector(462,11);
      WHEN "101100111" => data <= conv_std_logic_vector(461,11);
      WHEN "101101000" => data <= conv_std_logic_vector(460,11);
      WHEN "101101001" => data <= conv_std_logic_vector(460,11);
      WHEN "101101010" => data <= conv_std_logic_vector(459,11);
      WHEN "101101011" => data <= conv_std_logic_vector(458,11);
      WHEN "101101100" => data <= conv_std_logic_vector(457,11);
      WHEN "101101101" => data <= conv_std_logic_vector(456,11);
      WHEN "101101110" => data <= conv_std_logic_vector(456,11);
      WHEN "101101111" => data <= conv_std_logic_vector(455,11);
      WHEN "101110000" => data <= conv_std_logic_vector(454,11);
      WHEN "101110001" => data <= conv_std_logic_vector(453,11);
      WHEN "101110010" => data <= conv_std_logic_vector(453,11);
      WHEN "101110011" => data <= conv_std_logic_vector(452,11);
      WHEN "101110100" => data <= conv_std_logic_vector(451,11);
      WHEN "101110101" => data <= conv_std_logic_vector(450,11);
      WHEN "101110110" => data <= conv_std_logic_vector(449,11);
      WHEN "101110111" => data <= conv_std_logic_vector(449,11);
      WHEN "101111000" => data <= conv_std_logic_vector(448,11);
      WHEN "101111001" => data <= conv_std_logic_vector(447,11);
      WHEN "101111010" => data <= conv_std_logic_vector(446,11);
      WHEN "101111011" => data <= conv_std_logic_vector(446,11);
      WHEN "101111100" => data <= conv_std_logic_vector(445,11);
      WHEN "101111101" => data <= conv_std_logic_vector(444,11);
      WHEN "101111110" => data <= conv_std_logic_vector(443,11);
      WHEN "101111111" => data <= conv_std_logic_vector(443,11);
      WHEN "110000000" => data <= conv_std_logic_vector(442,11);
      WHEN "110000001" => data <= conv_std_logic_vector(441,11);
      WHEN "110000010" => data <= conv_std_logic_vector(440,11);
      WHEN "110000011" => data <= conv_std_logic_vector(440,11);
      WHEN "110000100" => data <= conv_std_logic_vector(439,11);
      WHEN "110000101" => data <= conv_std_logic_vector(438,11);
      WHEN "110000110" => data <= conv_std_logic_vector(438,11);
      WHEN "110000111" => data <= conv_std_logic_vector(437,11);
      WHEN "110001000" => data <= conv_std_logic_vector(436,11);
      WHEN "110001001" => data <= conv_std_logic_vector(435,11);
      WHEN "110001010" => data <= conv_std_logic_vector(435,11);
      WHEN "110001011" => data <= conv_std_logic_vector(434,11);
      WHEN "110001100" => data <= conv_std_logic_vector(433,11);
      WHEN "110001101" => data <= conv_std_logic_vector(433,11);
      WHEN "110001110" => data <= conv_std_logic_vector(432,11);
      WHEN "110001111" => data <= conv_std_logic_vector(431,11);
      WHEN "110010000" => data <= conv_std_logic_vector(430,11);
      WHEN "110010001" => data <= conv_std_logic_vector(430,11);
      WHEN "110010010" => data <= conv_std_logic_vector(429,11);
      WHEN "110010011" => data <= conv_std_logic_vector(428,11);
      WHEN "110010100" => data <= conv_std_logic_vector(428,11);
      WHEN "110010101" => data <= conv_std_logic_vector(427,11);
      WHEN "110010110" => data <= conv_std_logic_vector(426,11);
      WHEN "110010111" => data <= conv_std_logic_vector(425,11);
      WHEN "110011000" => data <= conv_std_logic_vector(425,11);
      WHEN "110011001" => data <= conv_std_logic_vector(424,11);
      WHEN "110011010" => data <= conv_std_logic_vector(423,11);
      WHEN "110011011" => data <= conv_std_logic_vector(423,11);
      WHEN "110011100" => data <= conv_std_logic_vector(422,11);
      WHEN "110011101" => data <= conv_std_logic_vector(421,11);
      WHEN "110011110" => data <= conv_std_logic_vector(421,11);
      WHEN "110011111" => data <= conv_std_logic_vector(420,11);
      WHEN "110100000" => data <= conv_std_logic_vector(419,11);
      WHEN "110100001" => data <= conv_std_logic_vector(419,11);
      WHEN "110100010" => data <= conv_std_logic_vector(418,11);
      WHEN "110100011" => data <= conv_std_logic_vector(417,11);
      WHEN "110100100" => data <= conv_std_logic_vector(417,11);
      WHEN "110100101" => data <= conv_std_logic_vector(416,11);
      WHEN "110100110" => data <= conv_std_logic_vector(415,11);
      WHEN "110100111" => data <= conv_std_logic_vector(415,11);
      WHEN "110101000" => data <= conv_std_logic_vector(414,11);
      WHEN "110101001" => data <= conv_std_logic_vector(413,11);
      WHEN "110101010" => data <= conv_std_logic_vector(413,11);
      WHEN "110101011" => data <= conv_std_logic_vector(412,11);
      WHEN "110101100" => data <= conv_std_logic_vector(411,11);
      WHEN "110101101" => data <= conv_std_logic_vector(411,11);
      WHEN "110101110" => data <= conv_std_logic_vector(410,11);
      WHEN "110101111" => data <= conv_std_logic_vector(409,11);
      WHEN "110110000" => data <= conv_std_logic_vector(409,11);
      WHEN "110110001" => data <= conv_std_logic_vector(408,11);
      WHEN "110110010" => data <= conv_std_logic_vector(407,11);
      WHEN "110110011" => data <= conv_std_logic_vector(407,11);
      WHEN "110110100" => data <= conv_std_logic_vector(406,11);
      WHEN "110110101" => data <= conv_std_logic_vector(405,11);
      WHEN "110110110" => data <= conv_std_logic_vector(405,11);
      WHEN "110110111" => data <= conv_std_logic_vector(404,11);
      WHEN "110111000" => data <= conv_std_logic_vector(404,11);
      WHEN "110111001" => data <= conv_std_logic_vector(403,11);
      WHEN "110111010" => data <= conv_std_logic_vector(402,11);
      WHEN "110111011" => data <= conv_std_logic_vector(402,11);
      WHEN "110111100" => data <= conv_std_logic_vector(401,11);
      WHEN "110111101" => data <= conv_std_logic_vector(400,11);
      WHEN "110111110" => data <= conv_std_logic_vector(400,11);
      WHEN "110111111" => data <= conv_std_logic_vector(399,11);
      WHEN "111000000" => data <= conv_std_logic_vector(399,11);
      WHEN "111000001" => data <= conv_std_logic_vector(398,11);
      WHEN "111000010" => data <= conv_std_logic_vector(397,11);
      WHEN "111000011" => data <= conv_std_logic_vector(397,11);
      WHEN "111000100" => data <= conv_std_logic_vector(396,11);
      WHEN "111000101" => data <= conv_std_logic_vector(395,11);
      WHEN "111000110" => data <= conv_std_logic_vector(395,11);
      WHEN "111000111" => data <= conv_std_logic_vector(394,11);
      WHEN "111001000" => data <= conv_std_logic_vector(394,11);
      WHEN "111001001" => data <= conv_std_logic_vector(393,11);
      WHEN "111001010" => data <= conv_std_logic_vector(392,11);
      WHEN "111001011" => data <= conv_std_logic_vector(392,11);
      WHEN "111001100" => data <= conv_std_logic_vector(391,11);
      WHEN "111001101" => data <= conv_std_logic_vector(391,11);
      WHEN "111001110" => data <= conv_std_logic_vector(390,11);
      WHEN "111001111" => data <= conv_std_logic_vector(389,11);
      WHEN "111010000" => data <= conv_std_logic_vector(389,11);
      WHEN "111010001" => data <= conv_std_logic_vector(388,11);
      WHEN "111010010" => data <= conv_std_logic_vector(388,11);
      WHEN "111010011" => data <= conv_std_logic_vector(387,11);
      WHEN "111010100" => data <= conv_std_logic_vector(386,11);
      WHEN "111010101" => data <= conv_std_logic_vector(386,11);
      WHEN "111010110" => data <= conv_std_logic_vector(385,11);
      WHEN "111010111" => data <= conv_std_logic_vector(385,11);
      WHEN "111011000" => data <= conv_std_logic_vector(384,11);
      WHEN "111011001" => data <= conv_std_logic_vector(383,11);
      WHEN "111011010" => data <= conv_std_logic_vector(383,11);
      WHEN "111011011" => data <= conv_std_logic_vector(382,11);
      WHEN "111011100" => data <= conv_std_logic_vector(382,11);
      WHEN "111011101" => data <= conv_std_logic_vector(381,11);
      WHEN "111011110" => data <= conv_std_logic_vector(381,11);
      WHEN "111011111" => data <= conv_std_logic_vector(380,11);
      WHEN "111100000" => data <= conv_std_logic_vector(379,11);
      WHEN "111100001" => data <= conv_std_logic_vector(379,11);
      WHEN "111100010" => data <= conv_std_logic_vector(378,11);
      WHEN "111100011" => data <= conv_std_logic_vector(378,11);
      WHEN "111100100" => data <= conv_std_logic_vector(377,11);
      WHEN "111100101" => data <= conv_std_logic_vector(377,11);
      WHEN "111100110" => data <= conv_std_logic_vector(376,11);
      WHEN "111100111" => data <= conv_std_logic_vector(375,11);
      WHEN "111101000" => data <= conv_std_logic_vector(375,11);
      WHEN "111101001" => data <= conv_std_logic_vector(374,11);
      WHEN "111101010" => data <= conv_std_logic_vector(374,11);
      WHEN "111101011" => data <= conv_std_logic_vector(373,11);
      WHEN "111101100" => data <= conv_std_logic_vector(373,11);
      WHEN "111101101" => data <= conv_std_logic_vector(372,11);
      WHEN "111101110" => data <= conv_std_logic_vector(372,11);
      WHEN "111101111" => data <= conv_std_logic_vector(371,11);
      WHEN "111110000" => data <= conv_std_logic_vector(370,11);
      WHEN "111110001" => data <= conv_std_logic_vector(370,11);
      WHEN "111110010" => data <= conv_std_logic_vector(369,11);
      WHEN "111110011" => data <= conv_std_logic_vector(369,11);
      WHEN "111110100" => data <= conv_std_logic_vector(368,11);
      WHEN "111110101" => data <= conv_std_logic_vector(368,11);
      WHEN "111110110" => data <= conv_std_logic_vector(367,11);
      WHEN "111110111" => data <= conv_std_logic_vector(367,11);
      WHEN "111111000" => data <= conv_std_logic_vector(366,11);
      WHEN "111111001" => data <= conv_std_logic_vector(366,11);
      WHEN "111111010" => data <= conv_std_logic_vector(365,11);
      WHEN "111111011" => data <= conv_std_logic_vector(364,11);
      WHEN "111111100" => data <= conv_std_logic_vector(364,11);
      WHEN "111111101" => data <= conv_std_logic_vector(363,11);
      WHEN "111111110" => data <= conv_std_logic_vector(363,11);
      WHEN "111111111" => data <= conv_std_logic_vector(362,11);
      WHEN others => data <= conv_std_logic_vector(0,11);
    END CASE;
  END PROCESS;

END rtl;