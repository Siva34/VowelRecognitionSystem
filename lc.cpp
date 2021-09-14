#define _USE_MATH_DEFINES
#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<math.h>
#include<vector>
#include <iomanip> 
using namespace std;

//dc shift calculation based on some initial silence samples and applied to entire input
void applyDCShift(vector<long double> &arr){
	long double dc=0;
	for(int i=0;i<1600;i++){
	   dc = dc+arr[i];
	}
	dc = dc/(double)1600;
	for(int i=0;i<arr.size();i++){
	   arr[i] = arr[i]+dc;
	}
}

//maximum allowed sample value is 10000 and accordingly every sample value is normalized based on normalized factor i.e., 10000/max_value
void performNormalization(vector<long double> &arr){
	long double max_val=-50000;
	for(int i=0;i<arr.size();i++){
	   if(arr[i] > max_val){
	      max_val = arr[i];
	   }
	}
	long double norm_fact;
	norm_fact = 10000/max_val;
	for(int i=0;i<arr.size();i++){
          arr[i] = arr[i]*norm_fact;		
	}
}

//once steady frames are noted down, hamming window is applied to every frame
void applyHammingWindow(vector<long double> &arr){
	int frames = arr.size()/320;
	long double h_w[320];
	for(int i=0;i<320;i++){
	    h_w[i] = 0.54 -( 0.46*(cos((2*M_PI*i)/double(319))));
	}
	for(int i = 0; i<frames ; i++){
	   long double val = 0;
	   int start = 320*i;
	   for(int j = 320*i ;j<(320*(i+1));j++){
	      arr[j] = arr[j]*h_w[j-start];
	   }
	 }
}

//this function returns maximum ste frame number and also it calculated every frame ste and store these values into ste[]
int calculateSTE(vector<long double> &arr,vector<long double> &ste){
	int frames = arr.size()/320;
	long double max_ste_index = 0;
	for(int i = 0; i<frames ; i++){
	   long double val = 0;
	   for(int j = 320*i ;j<(320*(i+1));j++){
	      val = val + arr[j]*arr[j];
	   }
	   val=val/320;
	   ste.push_back(val);
	   if(ste[i] >= ste[max_ste_index]){
	      max_ste_index = i;
	   }
	}
	return max_ste_index;
}

//ceps are calculated for the passed frame number by durbins method
void findCepsForThisFrame(long double* cep,long double* ai,long double* ri,vector<long double> &arr,int frame){
   	//Auto-correlation analysis
	int p = 12;
	int start = frame*320;  // starting index of the passed frame in arr
	int end = start+319;    // ending index of the passed frame in arr
	long double* r = new long double[13];
	for(int k=0;k<13;k++){
		    long double val = 0;
	        for(int m = start;m<=(end-k);m++){
			   val =  val + (arr[m]*arr[m+k]);
			}
			r[k] = val/320;
	}
	//form toeplitz matrix
    long double toe[12][12];
	for(int i=0;i<12;i++){
	     for(int j=i;j>=0;j--){
		     toe[i][j] = r[j];
		 }
		 for(int k=1;k<(12-i);k++){
		     toe[i][k] = r[k];
		 }
	}

	//durbins method
	long double E[13],k[13],a[13][13];
	E[0]=r[0];
	for(int i=1;i<=12;i++){
	   long double temp=0;
	   for(int j=1;j<=(i-1);j++){
	       temp = temp + (a[j][i-1]*r[i-j]);
	   }
	   k[i] = (r[i]-temp)/E[i-1];
	   a[i][i]=k[i];
	   for(int j=1;j<=(i-1);j++){
	       a[j][i] = a[j][i-1]-(k[i]*a[i-j][i-1]);
	   }
	   E[i]=(1-(k[i]*k[i]))*E[i-1];
	}
	long double* coeff = new long double[13];
	coeff[0] = 0;
	for(int i=1;i<=12;i++){
	   coeff[i] = a[i][12];
	}
    long double cepstral[13];
	cepstral[0] = log(r[0]*r[0]);
	for(int m=1;m<=12;m++){
	   long double temp = 0;
	   for(int k=1;k<=(m-1);k++){
	       temp = temp + ((k/(double)m)*cepstral[k]*(coeff[m-k]));
	   }
	    cepstral[m] = coeff[m] + temp;
	}
	for(int i=0;i<=12;i++){
	   ai[i] = coeff[i];
	   ri[i] = r[0];
	   cep[i] = cepstral[i];
	}
}

//applying raised sine window once ceps are calculated 
void applyRaisedSineWindow(long double* cep){
	int q = 12;
	for(int i=1;i<=12;i++){
	   long double w = 1 + ( (q/2)*sin((M_PI*i)/(double)q));
	   cep[i] = cep[i]*w;
	}
}

/*
trainFunc() helps in calculating cepstral coefficients of 5 steady frames for each vowel in following manner : 
for each vowel
1) for all 10 files of each vowel, we are going to calculate 5 steady frames ceps. i.e., in total 50 rows of ceps (in which 10 rows of ceps for
   each steady frame among 5)
2) calculate average values of ceps for 10 rows of each steady frame
3) now we will get 5 rows of ceps for 5 steady frames
these 5 rows of ceps are going to store into textfile named as vowel.txt (example = a.txt)
*/
void trainFunc(long double*** cep)
{
	string vowels[5] = {"a","e","i","o","u"};
	for(int v=0;v<5;v++){
		long double** temp = new long double*[51];  // temp is used for storing all 10*5 = 50 rows of steady frames ceps

	    for(int vow_rep=1;vow_rep<=10;vow_rep++){
	    	int start_frame_temp = (vow_rep-1)*5+1;
	    	int end_frame_temp = start_frame_temp+4;

	    	for(int i=start_frame_temp;i<=end_frame_temp;i++){
	    		temp[i] = new long double[13];     // as each row in temp is having 12 columns of ceps from c1 to c12
			}
			stringstream ss;
		    int roll_no = 204101037;
	        ss << roll_no;
		    string filename = ss.str()+"_";
		    filename += vowels[v];
		    filename += "_";
		    ss.str("");
		    ss.clear();
		    ss << vow_rep;
		    filename += ss.str();
		    filename += ".txt";
		    
		    fstream fs;
		    fs.open(filename,ios::in);
		    vector<long double> arr;
		    string sample;
		    while(getline(fs,sample)){
		    	arr.push_back(stold(sample));
			}
			fs.close();
			applyDCShift(arr);                            // DC SHIFT
		    performNormalization(arr);                    // NORMALIZATION
	        int frames = arr.size()/320;
	 	    vector<long double> ste;
		    int max_ste_frame = calculateSTE(arr,ste);    // MAX STE FRAME for finding steady frames
			applyHammingWindow(arr);                      // hamming window

		    int start_frame = max_ste_frame-2;
		    int end_frame = max_ste_frame+2;
			long double** ceps_steady_frames = new long double*[5];
		    long double** ai = new long double*[5];
		    long double** ri = new long double*[5];
		    for(int k = start_frame;k<=end_frame;k++){
		          ceps_steady_frames[k-start_frame] = new long double[13];
		          ai[k-start_frame] = new long double[13];
		          ri[k-start_frame] = new long double[13];
		          findCepsForThisFrame(ceps_steady_frames[k-start_frame],ai[k-start_frame],ri[k-start_frame],arr,k-start_frame);
		          applyRaisedSineWindow(ceps_steady_frames[k-start_frame]);   // applying raised sine window once ceps are calculated
		   }
		   
		   	for(int i=0;i<5;i++){
			    for(int j=1;j<=12;j++){
				    temp[start_frame_temp+i][j] = ceps_steady_frames[i][j];
			    }
			} 
		   
		}
		for(int frame=1;frame<=5;frame++){
			long double c[13]={0};
			for(int i=frame;i<=50;i=i+5){
				for(int k=1;k<=12;k++){
					c[k] = c[k] + temp[i][k];
				}
			}
			for(int i=1;i<=12;i++){
				cep[v][frame-1][i] = c[i]/10;
			}
		}
		
	}
	
	for(int i=0;i<5;i++){
	   cout<<"for vowel : "<<vowels[i]<<endl;
	   fstream t;
	   string outputfile;
	   outputfile = vowels[i] +".txt";
	   t.open(outputfile,ios::out);
	   for(int j=0;j<5;j++){
		   cout<<"for frame "<<j+1<<" ceps are : "<<endl;
	       for(int k=1;k<=12;k++){
		       t << cep[i][j][k]<<endl; // storing cepstral coefficients of all 5 frames for a vowel into text file line by line
			   cout<<cep[i][j][k] <<" "; 
		   }
		   cout<<endl;
	   }
	   cout<<endl;
	}

}

//this function returns which vowel we speak out based on tokhura distance between trained ceps and the tested data ceps going to calculate in this method
string testBasedOnCeps(vector<long double> &arr_test,string* vow){

	    applyDCShift(arr_test);
		performNormalization(arr_test);
		vector<long double> ste;
		int max_ste_frame = calculateSTE(arr_test,ste);
		applyHammingWindow(arr_test);
		int start_frame = max_ste_frame-2;
		int end_frame = max_ste_frame+2;
		long double** ceps_steady_frames = new long double*[5];
		long double** ai = new long double*[5];
		long double** ri = new long double*[5];
		for(int k = start_frame;k<=end_frame;k++){
		   ceps_steady_frames[k-start_frame] = new long double[13];
		   ai[k-start_frame] = new long double[13];
		   ri[k-start_frame] = new long double[13];
		   findCepsForThisFrame(ceps_steady_frames[k-start_frame],ai[k-start_frame],ri[k-start_frame],arr_test,k-start_frame);
		   applyRaisedSineWindow(ceps_steady_frames[k-start_frame]);
		}
	
	    vector<long double> tokhuraDist;
		vector<double> weights;
		weights.push_back(1.0);     // these weights are provided by sir
		weights.push_back(3.0);
		weights.push_back(7.0);
		weights.push_back(13.0);
		weights.push_back(19.0);
		weights.push_back(22.0);
		weights.push_back(25.0);
		weights.push_back(33.0);
		weights.push_back(42.0);
		weights.push_back(50.0);
		weights.push_back(56.0);
        weights.push_back(61.0);
		
		for(int i=0;i<5;i++){
		   string str = vow[i]+".txt";
		   vector<long double> temp;
		   fstream inp;
		   string cep_v;
		   inp.open(str,ios::in);
		   while(getline(inp,cep_v)){
			   temp.push_back(stold(cep_v));
		   }
		   inp.close();

		   long double **cep_r = new long double*[5];
		   for(int j=0;j<5;j++){
		       cep_r[j] = new long double[13];
			   for(int k=1;k<=12;k++){
			      cep_r[j][k] = temp[(j*12)+(k-1)];
			   }

		   }
		   long double dist_frames[5]={0};
		   //calculating distances of all 5 frames and their sum is going to store into tokhuraDist vector array
		   for(int frame=0;frame<5;frame++){
			    long double dist_val = 0;
		        for(int i=1;i<=12;i++){
		            dist_val = dist_val + ((weights[i-1])*((cep_r[frame][i]-ceps_steady_frames[frame][i])*(cep_r[frame][i]-ceps_steady_frames[frame][i])));
		        }
				dist_frames[frame] = dist_val;
		   }
		   long double x=0;
		   for(int i=0;i<5;i++){
		      x += dist_frames[i];
		   }
		   x=x/5;
		   tokhuraDist.push_back(x);
		   cout<<"tokhura distance = "<<x<<endl;
		   temp.clear();

		 }
		  // once all 5 vowels to testing_vowel tokhura distances are calculated, decide vowel spoken based on minimum distance among 5
		  long double min_dist = 123445525;
		  int index = -1;
		  for(int i=0;i<5;i++){
		     if(tokhuraDist[i]<min_dist){
		          min_dist = tokhuraDist[i];
				  index = i;
		     }
		  }
		 cout<<endl<<"recognized one = "<<vow[index]<<endl;
		 tokhuraDist.clear();
		 weights.clear();
		 arr_test.clear();
		 return vow[index];   // returning which vowel we recognized 
}

//testFunc() used for testing pre-recorded text files(11 to 20 of each file) and among 50 files passed, how many are recognizing correctly is calculated
void testFunc(long double*** cep){
	string vow[5] = {"a","e","i","o","u"};
	int count = 0;
	for(int i=0;i<5;i++){
		for(int j=11;j<=20;j++){
        string str = "204101037_"+vow[i];
		str = str+"_";
	    stringstream ss;
	    ss << j;
		str += ss.str();
		str += ".txt";
		ss.str("");
		ss.clear();
		fstream f;
		cout<<endl<<str<<endl;
		f.open(str,ios::in);
		string sample;
		vector<long double> arr_test;
		while(getline(f, sample)){
			arr_test.push_back(stold(sample));
	    }
	    f.close();
		string res = testBasedOnCeps(arr_test,vow);
		if(res==vow[i]){
		   count++;
		}
		
		}
	}
	cout<<endl<<"out of 50, "<<count<<" are passed"<<endl;
}

//this is used for vowel recognition through live-recording module instead of pre-recorded text files	
void testWithYourVoice(){
   //record the vowel
	string vow[5] = {"a","e","i","o","u"};
	bool flag = true;
	while(flag){
	       system("Recording_Module.exe 5 input_file.wav voice.txt"); 
           string str = "voice.txt";
	       fstream t;
	       t.open(str,ios::in);
	       vector<long double> arr_t;
	       string ip;
	       int count=0;
	       while(getline(t,ip)){
	              count++;
	              if(count<16000 || count>=48000){
	                   continue;
	              }
	              arr_t.push_back(stold(ip));
	       }
	      t.close();
	      string ouput = testBasedOnCeps(arr_t,vow);
		  int num;
		  cout << "do you want to test more ? "<<endl;
		  cout << "choose 1 for yes and 2 for no = ";
		  cin >> num;
		  if(num == 2){
		      flag = false;
		  }

   }
}

int main(){
	cout<<"training function for vowel recoginition is started "<<endl;
	long double*** cep = new long double**[5];
	for(int i=0;i<5;i++){
	   cep[i] = new long double*[5];
	   for(int j=0;j<5;j++){
	   	  cep[i][j] = new long double[13];
	   }
	}
	trainFunc(cep);
	cout<<endl;
	cout<<"how you want to perform testing = ";
	cout<<"1) pre-recorded 2) live-recording"<<endl;
	cout<<"please enter your choice = ";
	int choice;
	cin >> choice;
	if(choice == 1){
	     testFunc(cep);
	}
	else{
	     testWithYourVoice();
	}
	return 0;
}

