# VowelRecognitionSystem
========================================================================
    CONSOLE APPLICATION : lc Project Overview
========================================================================

cep is a 3 dimensional array, and can interpret like below:
         where cep[0] refers to vowel a
               cep[1] refers to vowel e
 
           and cep[0][0] points to 1st steady frame of vowel a
               cep[2][3] points to 4th steady frame of vowel i
	 
           and cep[0][0][1] points to c1 of 1st steady frame of vowel a
               cep[2][3][12] points to c12 of 4th steady frame of vowel i

here i am doing 2 things, 
step 1 : training the vowel recognition system in following manner :
   
         trainFunc() helps in calculating cepstral coefficients of 5 steady frames for each vowel in following manner : 
         for each vowel
         1) for all 10 files of each vowel, we are going to calculate 5 steady frames ceps. i.e., in total 50 rows of ceps (in which 10 rows of ceps for
            each steady frame among 5)
         2) calculate average values of ceps for 10 rows of each steady frame
         3) now we will get 5 rows of ceps for 5 steady frames
            these 5 rows of ceps are going to store into textfile named as vowel.txt (example = a.txt) 

step 2: testing the vowel recognition system

        here i implemented it through 2 modes, one is based on pre-recorded data and another through live-recording module 

	 	  	  
procedure for finding ceps of testing data through live recording :-
step 1 :- record one vowel utterence in 5 seconds time span
step 2 :- removing 1 second from both beginning and ending
step 3 :- apply dcshift
step 4 :- perform Normalization
step 5 :- finding maximum ste frame and note down steady frames
step 6 :- applying hamming window for each frame
step 7 :- by using auto-correlation analysis and durbins method, find out cepstral coefficients
step 8 :- apply raised sine window for output of step 7

once ceps are calculated, by finding out tokhura distances of each vowel to this recorded data, we are going to recognize vowel spoken (which is 
having minimum tokhura distance)
