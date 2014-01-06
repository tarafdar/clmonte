
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
--***   FP_EXPLUT7.VHD                            ***
--***                                             ***
--***   Function: Look Up Table - EXP()           ***
--***                                             ***
--***   Generated by MATLAB Utility               ***
--***                                             ***
--***   18/02/08 ML                               ***
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

ENTITY fp_explut7 IS
PORT (
      address : IN STD_LOGIC_VECTOR (7 DOWNTO 1);
      mantissa : OUT STD_LOGIC_VECTOR (23 DOWNTO 1);
      exponent : OUT STD_LOGIC
     );
END fp_explut7;

ARCHITECTURE rtl OF fp_explut7 IS

BEGIN

  pca: PROCESS (address)
  BEGIN
    CASE address IS
      WHEN "0000000" =>
            mantissa <= conv_std_logic_vector(0,23);
            exponent <= '0';
      WHEN "0000001" =>
            mantissa <= conv_std_logic_vector(65793,23);
            exponent <= '0';
      WHEN "0000010" =>
            mantissa <= conv_std_logic_vector(132101,23);
            exponent <= '0';
      WHEN "0000011" =>
            mantissa <= conv_std_logic_vector(198930,23);
            exponent <= '0';
      WHEN "0000100" =>
            mantissa <= conv_std_logic_vector(266283,23);
            exponent <= '0';
      WHEN "0000101" =>
            mantissa <= conv_std_logic_vector(334164,23);
            exponent <= '0';
      WHEN "0000110" =>
            mantissa <= conv_std_logic_vector(402578,23);
            exponent <= '0';
      WHEN "0000111" =>
            mantissa <= conv_std_logic_vector(471528,23);
            exponent <= '0';
      WHEN "0001000" =>
            mantissa <= conv_std_logic_vector(541019,23);
            exponent <= '0';
      WHEN "0001001" =>
            mantissa <= conv_std_logic_vector(611055,23);
            exponent <= '0';
      WHEN "0001010" =>
            mantissa <= conv_std_logic_vector(681640,23);
            exponent <= '0';
      WHEN "0001011" =>
            mantissa <= conv_std_logic_vector(752779,23);
            exponent <= '0';
      WHEN "0001100" =>
            mantissa <= conv_std_logic_vector(824476,23);
            exponent <= '0';
      WHEN "0001101" =>
            mantissa <= conv_std_logic_vector(896735,23);
            exponent <= '0';
      WHEN "0001110" =>
            mantissa <= conv_std_logic_vector(969560,23);
            exponent <= '0';
      WHEN "0001111" =>
            mantissa <= conv_std_logic_vector(1042957,23);
            exponent <= '0';
      WHEN "0010000" =>
            mantissa <= conv_std_logic_vector(1116930,23);
            exponent <= '0';
      WHEN "0010001" =>
            mantissa <= conv_std_logic_vector(1191483,23);
            exponent <= '0';
      WHEN "0010010" =>
            mantissa <= conv_std_logic_vector(1266621,23);
            exponent <= '0';
      WHEN "0010011" =>
            mantissa <= conv_std_logic_vector(1342348,23);
            exponent <= '0';
      WHEN "0010100" =>
            mantissa <= conv_std_logic_vector(1418668,23);
            exponent <= '0';
      WHEN "0010101" =>
            mantissa <= conv_std_logic_vector(1495588,23);
            exponent <= '0';
      WHEN "0010110" =>
            mantissa <= conv_std_logic_vector(1573110,23);
            exponent <= '0';
      WHEN "0010111" =>
            mantissa <= conv_std_logic_vector(1651241,23);
            exponent <= '0';
      WHEN "0011000" =>
            mantissa <= conv_std_logic_vector(1729985,23);
            exponent <= '0';
      WHEN "0011001" =>
            mantissa <= conv_std_logic_vector(1809346,23);
            exponent <= '0';
      WHEN "0011010" =>
            mantissa <= conv_std_logic_vector(1889329,23);
            exponent <= '0';
      WHEN "0011011" =>
            mantissa <= conv_std_logic_vector(1969940,23);
            exponent <= '0';
      WHEN "0011100" =>
            mantissa <= conv_std_logic_vector(2051183,23);
            exponent <= '0';
      WHEN "0011101" =>
            mantissa <= conv_std_logic_vector(2133064,23);
            exponent <= '0';
      WHEN "0011110" =>
            mantissa <= conv_std_logic_vector(2215586,23);
            exponent <= '0';
      WHEN "0011111" =>
            mantissa <= conv_std_logic_vector(2298756,23);
            exponent <= '0';
      WHEN "0100000" =>
            mantissa <= conv_std_logic_vector(2382578,23);
            exponent <= '0';
      WHEN "0100001" =>
            mantissa <= conv_std_logic_vector(2467057,23);
            exponent <= '0';
      WHEN "0100010" =>
            mantissa <= conv_std_logic_vector(2552199,23);
            exponent <= '0';
      WHEN "0100011" =>
            mantissa <= conv_std_logic_vector(2638009,23);
            exponent <= '0';
      WHEN "0100100" =>
            mantissa <= conv_std_logic_vector(2724492,23);
            exponent <= '0';
      WHEN "0100101" =>
            mantissa <= conv_std_logic_vector(2811653,23);
            exponent <= '0';
      WHEN "0100110" =>
            mantissa <= conv_std_logic_vector(2899498,23);
            exponent <= '0';
      WHEN "0100111" =>
            mantissa <= conv_std_logic_vector(2988032,23);
            exponent <= '0';
      WHEN "0101000" =>
            mantissa <= conv_std_logic_vector(3077260,23);
            exponent <= '0';
      WHEN "0101001" =>
            mantissa <= conv_std_logic_vector(3167188,23);
            exponent <= '0';
      WHEN "0101010" =>
            mantissa <= conv_std_logic_vector(3257821,23);
            exponent <= '0';
      WHEN "0101011" =>
            mantissa <= conv_std_logic_vector(3349165,23);
            exponent <= '0';
      WHEN "0101100" =>
            mantissa <= conv_std_logic_vector(3441225,23);
            exponent <= '0';
      WHEN "0101101" =>
            mantissa <= conv_std_logic_vector(3534008,23);
            exponent <= '0';
      WHEN "0101110" =>
            mantissa <= conv_std_logic_vector(3627518,23);
            exponent <= '0';
      WHEN "0101111" =>
            mantissa <= conv_std_logic_vector(3721762,23);
            exponent <= '0';
      WHEN "0110000" =>
            mantissa <= conv_std_logic_vector(3816745,23);
            exponent <= '0';
      WHEN "0110001" =>
            mantissa <= conv_std_logic_vector(3912472,23);
            exponent <= '0';
      WHEN "0110010" =>
            mantissa <= conv_std_logic_vector(4008951,23);
            exponent <= '0';
      WHEN "0110011" =>
            mantissa <= conv_std_logic_vector(4106186,23);
            exponent <= '0';
      WHEN "0110100" =>
            mantissa <= conv_std_logic_vector(4204184,23);
            exponent <= '0';
      WHEN "0110101" =>
            mantissa <= conv_std_logic_vector(4302951,23);
            exponent <= '0';
      WHEN "0110110" =>
            mantissa <= conv_std_logic_vector(4402492,23);
            exponent <= '0';
      WHEN "0110111" =>
            mantissa <= conv_std_logic_vector(4502814,23);
            exponent <= '0';
      WHEN "0111000" =>
            mantissa <= conv_std_logic_vector(4603922,23);
            exponent <= '0';
      WHEN "0111001" =>
            mantissa <= conv_std_logic_vector(4705824,23);
            exponent <= '0';
      WHEN "0111010" =>
            mantissa <= conv_std_logic_vector(4808525,23);
            exponent <= '0';
      WHEN "0111011" =>
            mantissa <= conv_std_logic_vector(4912031,23);
            exponent <= '0';
      WHEN "0111100" =>
            mantissa <= conv_std_logic_vector(5016349,23);
            exponent <= '0';
      WHEN "0111101" =>
            mantissa <= conv_std_logic_vector(5121486,23);
            exponent <= '0';
      WHEN "0111110" =>
            mantissa <= conv_std_logic_vector(5227447,23);
            exponent <= '0';
      WHEN "0111111" =>
            mantissa <= conv_std_logic_vector(5334239,23);
            exponent <= '0';
      WHEN "1000000" =>
            mantissa <= conv_std_logic_vector(5441868,23);
            exponent <= '0';
      WHEN "1000001" =>
            mantissa <= conv_std_logic_vector(5550342,23);
            exponent <= '0';
      WHEN "1000010" =>
            mantissa <= conv_std_logic_vector(5659667,23);
            exponent <= '0';
      WHEN "1000011" =>
            mantissa <= conv_std_logic_vector(5769849,23);
            exponent <= '0';
      WHEN "1000100" =>
            mantissa <= conv_std_logic_vector(5880895,23);
            exponent <= '0';
      WHEN "1000101" =>
            mantissa <= conv_std_logic_vector(5992812,23);
            exponent <= '0';
      WHEN "1000110" =>
            mantissa <= conv_std_logic_vector(6105607,23);
            exponent <= '0';
      WHEN "1000111" =>
            mantissa <= conv_std_logic_vector(6219286,23);
            exponent <= '0';
      WHEN "1001000" =>
            mantissa <= conv_std_logic_vector(6333858,23);
            exponent <= '0';
      WHEN "1001001" =>
            mantissa <= conv_std_logic_vector(6449327,23);
            exponent <= '0';
      WHEN "1001010" =>
            mantissa <= conv_std_logic_vector(6565703,23);
            exponent <= '0';
      WHEN "1001011" =>
            mantissa <= conv_std_logic_vector(6682991,23);
            exponent <= '0';
      WHEN "1001100" =>
            mantissa <= conv_std_logic_vector(6801199,23);
            exponent <= '0';
      WHEN "1001101" =>
            mantissa <= conv_std_logic_vector(6920334,23);
            exponent <= '0';
      WHEN "1001110" =>
            mantissa <= conv_std_logic_vector(7040403,23);
            exponent <= '0';
      WHEN "1001111" =>
            mantissa <= conv_std_logic_vector(7161415,23);
            exponent <= '0';
      WHEN "1010000" =>
            mantissa <= conv_std_logic_vector(7283375,23);
            exponent <= '0';
      WHEN "1010001" =>
            mantissa <= conv_std_logic_vector(7406292,23);
            exponent <= '0';
      WHEN "1010010" =>
            mantissa <= conv_std_logic_vector(7530173,23);
            exponent <= '0';
      WHEN "1010011" =>
            mantissa <= conv_std_logic_vector(7655025,23);
            exponent <= '0';
      WHEN "1010100" =>
            mantissa <= conv_std_logic_vector(7780857,23);
            exponent <= '0';
      WHEN "1010101" =>
            mantissa <= conv_std_logic_vector(7907676,23);
            exponent <= '0';
      WHEN "1010110" =>
            mantissa <= conv_std_logic_vector(8035489,23);
            exponent <= '0';
      WHEN "1010111" =>
            mantissa <= conv_std_logic_vector(8164305,23);
            exponent <= '0';
      WHEN "1011000" =>
            mantissa <= conv_std_logic_vector(8294131,23);
            exponent <= '0';
      WHEN "1011001" =>
            mantissa <= conv_std_logic_vector(18184,23);
            exponent <= '1';
      WHEN "1011010" =>
            mantissa <= conv_std_logic_vector(84119,23);
            exponent <= '1';
      WHEN "1011011" =>
            mantissa <= conv_std_logic_vector(150571,23);
            exponent <= '1';
      WHEN "1011100" =>
            mantissa <= conv_std_logic_vector(217545,23);
            exponent <= '1';
      WHEN "1011101" =>
            mantissa <= conv_std_logic_vector(285044,23);
            exponent <= '1';
      WHEN "1011110" =>
            mantissa <= conv_std_logic_vector(353072,23);
            exponent <= '1';
      WHEN "1011111" =>
            mantissa <= conv_std_logic_vector(421634,23);
            exponent <= '1';
      WHEN "1100000" =>
            mantissa <= conv_std_logic_vector(490734,23);
            exponent <= '1';
      WHEN "1100001" =>
            mantissa <= conv_std_logic_vector(560375,23);
            exponent <= '1';
      WHEN "1100010" =>
            mantissa <= conv_std_logic_vector(630563,23);
            exponent <= '1';
      WHEN "1100011" =>
            mantissa <= conv_std_logic_vector(701301,23);
            exponent <= '1';
      WHEN "1100100" =>
            mantissa <= conv_std_logic_vector(772594,23);
            exponent <= '1';
      WHEN "1100101" =>
            mantissa <= conv_std_logic_vector(844446,23);
            exponent <= '1';
      WHEN "1100110" =>
            mantissa <= conv_std_logic_vector(916862,23);
            exponent <= '1';
      WHEN "1100111" =>
            mantissa <= conv_std_logic_vector(989846,23);
            exponent <= '1';
      WHEN "1101000" =>
            mantissa <= conv_std_logic_vector(1063402,23);
            exponent <= '1';
      WHEN "1101001" =>
            mantissa <= conv_std_logic_vector(1137535,23);
            exponent <= '1';
      WHEN "1101010" =>
            mantissa <= conv_std_logic_vector(1212249,23);
            exponent <= '1';
      WHEN "1101011" =>
            mantissa <= conv_std_logic_vector(1287550,23);
            exponent <= '1';
      WHEN "1101100" =>
            mantissa <= conv_std_logic_vector(1363441,23);
            exponent <= '1';
      WHEN "1101101" =>
            mantissa <= conv_std_logic_vector(1439927,23);
            exponent <= '1';
      WHEN "1101110" =>
            mantissa <= conv_std_logic_vector(1517013,23);
            exponent <= '1';
      WHEN "1101111" =>
            mantissa <= conv_std_logic_vector(1594704,23);
            exponent <= '1';
      WHEN "1110000" =>
            mantissa <= conv_std_logic_vector(1673004,23);
            exponent <= '1';
      WHEN "1110001" =>
            mantissa <= conv_std_logic_vector(1751918,23);
            exponent <= '1';
      WHEN "1110010" =>
            mantissa <= conv_std_logic_vector(1831452,23);
            exponent <= '1';
      WHEN "1110011" =>
            mantissa <= conv_std_logic_vector(1911608,23);
            exponent <= '1';
      WHEN "1110100" =>
            mantissa <= conv_std_logic_vector(1992394,23);
            exponent <= '1';
      WHEN "1110101" =>
            mantissa <= conv_std_logic_vector(2073813,23);
            exponent <= '1';
      WHEN "1110110" =>
            mantissa <= conv_std_logic_vector(2155871,23);
            exponent <= '1';
      WHEN "1110111" =>
            mantissa <= conv_std_logic_vector(2238572,23);
            exponent <= '1';
      WHEN "1111000" =>
            mantissa <= conv_std_logic_vector(2321922,23);
            exponent <= '1';
      WHEN "1111001" =>
            mantissa <= conv_std_logic_vector(2405926,23);
            exponent <= '1';
      WHEN "1111010" =>
            mantissa <= conv_std_logic_vector(2490589,23);
            exponent <= '1';
      WHEN "1111011" =>
            mantissa <= conv_std_logic_vector(2575915,23);
            exponent <= '1';
      WHEN "1111100" =>
            mantissa <= conv_std_logic_vector(2661911,23);
            exponent <= '1';
      WHEN "1111101" =>
            mantissa <= conv_std_logic_vector(2748582,23);
            exponent <= '1';
      WHEN "1111110" =>
            mantissa <= conv_std_logic_vector(2835932,23);
            exponent <= '1';
      WHEN "1111111" =>
            mantissa <= conv_std_logic_vector(2923967,23);
            exponent <= '1';
      WHEN others =>
           mantissa <= conv_std_logic_vector(0,23);
           exponent <= '0';
    END CASE;
  END PROCESS;

END rtl;
