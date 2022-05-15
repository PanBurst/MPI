
#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <chrono>
#include <vector>
#include <complex>
#include <mpi.h>

using std::cout;
using std::string;
using std::complex;
using std::vector;


#define TAG_RESULT 0

#define TAG_INIT 1
#define TAG_UPDATE 2
#define TAG_JOB_FINISHED 3
#define TAG_NEW_JOB 4
#define TAG_NO_NEW_JOB 5
#define TAG_TERMINATE 6

#define updateWsize 16



double calcMandelFracValue(double x, double y, int maxIters) {
	int i = 1;
	double reColor;
	double r, g, b, a;
	r = g = b = 0;
	a = 255;
	complex<double> _c = complex<double>(x, y);
	complex<double> d = complex<double>(2.0, 0.0);
	complex<double> z1 = complex<double>(0.0, 0.0);
	for (; i <= maxIters; i += 1) {
		//z1 = std::pow(z1, d) + _c;
		z1 = z1 * z1 + _c;
		if (std::abs(z1) > 2) {

			break;
		}
	}
	if (i < maxIters) {
		double v = std::log(i + 1.5 + std::log2(std::log(std::abs(z1)))) / 3.4;
		if (v < 1.0) {
			r = v * v * v * v;
			g = pow(v, 2.5);
			b = v;
		}
		else {
			r = v;
			g = pow(v, 1.5);
			b = v * v * v;
		}
	}
	r *= 255;
	g *= 255;
	b *= 255;
	reColor = round(r);
	reColor = round(g);
	reColor = round(b);
	reColor = round(a);
	return reColor;
}

int argc_;
char** argv_;

class RangeParameters {
public:
	int x, y;
	RangeParameters(int x, int y) {
		this->x = x;
		this->y = y;
	}
	int GetWork() {
		return y - x;
	}
};
class ControllerEqualRange {
	int w, h, cores;
	vector<RangeParameters*> parameters;
public:
	ControllerEqualRange(int h, int w, int cores) {
		this->h = h;
		this->w = w;
		this->cores = cores;
		int wStep = w / cores;
		parameters = vector<RangeParameters*>();
		for (int core = 0; core < cores; core += 1) {
			RangeParameters* param = new RangeParameters(core * wStep, (core + 1) * wStep);
			parameters.push_back(param);
		}
	}
	~ControllerEqualRange() {
		for (int core = 0; core < cores; core += 1) {
			delete parameters[cores];
		}
	}
	void Update(int x, int y, int core) {
		parameters[core - 1]->x = x;
		parameters[core - 1]->y = y;
	}
	RangeParameters* GetWork(int core) {
		if (parameters[core - 1]->GetWork() > 0) {
			return parameters[core - 1];
		}
		for (int id = 0; id < cores; id++) {
			int workLeft = parameters[id]->GetWork();
			if (workLeft >= updateWsize*2) {
				workLeft /= 2;
				parameters[core - 1]->x = parameters[id]->y - workLeft;
				parameters[core - 1]->y = parameters[id]->y;
				parameters[id]->y -= workLeft;
				return parameters[core - 1];
			}
		}
		return NULL;
	}
	bool Finished() {
		for (int id = 0; id < cores; id++) {
			if (parameters[id]->GetWork() >> 0) {
				return false;
			}
		}
		return true;
	}
	void print() {
		for (int i = 0; i < parameters.size(); i++) {
			cout << "Core: " << i + 1 << std::endl;
			cout << "\tRange: " << parameters[i]->x << "-" << parameters[i]->y << std::endl;
		}

	}


};





class Window {

public:
	int h, w;
	vector<int> sourceWaiting;
	int** pixels;
	Window(int h, int w) {
		this->h = h;
		this->w = w;
		pixels = new int* [w];
		for (int i = 0; i < w; i++) {
			pixels[i] = new int[h];
		}
		sourceWaiting = vector<int>();
	}
	~Window() {
		for (int i = 0; i < w; i++) {
			delete[] pixels[i];
		}
		delete[] pixels;
	}
	void changeRes(int h, int w) {
		this->h = h;
		this->w = w;
	}
	void close() {
		MPI_Status stat, stat2;
		for (auto& coreID : sourceWaiting) {
			double* buff = NULL;
			MPI_Send(buff, 0, MPI_DOUBLE, coreID, TAG_TERMINATE, MPI_COMM_WORLD);
		}
	}
	void RunTest() {
		double xr1 = -1.0, xr2 = 1.0, yr1 = -1.0, yr2 = 1.0;
		double xRange = abs(xr1) + abs(xr2), yRange = abs(yr1) + abs(yr2);
		double zoom = 2.0;
		//Event handler
		double realX = 0.0;
		double realY = 0.0;

		int samples = 10;
		double targetX = -0.45;
		double targetY = -0.8;
		int maxIters = 200;
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		//cout << "Frame,MaxIterations,SingleThread,EqualThreadRange,DynamicThread\n";
		cout << "RANK: " << rank << std::endl;
		if (rank == 0)
			cout << "Frame,MaxIterations,EqualThreadRange\n";


		for (int i = 1; i <= samples; i++) {
			double xStep = abs(targetX - realX) / (2.0 * i);
			double yStep = abs(targetY - realY) / (zoom * i);
			xStep = targetX > realX ? xStep : -xStep;
			yStep = targetY > realY ? yStep : -yStep;
			realX += xStep;
			realY += yStep;

			double nXRange = xRange / zoom;
			double nYRange = yRange / zoom;
			xr1 = realX - nXRange;
			xr2 = realX + nXRange;
			yr1 = realY - nYRange;
			yr2 = realY + nYRange;
			xRange = nXRange;
			yRange = nYRange;

			for (int iterations = 20; iterations <= maxIters; iterations += 20) {
				if (rank == 0)
					cout << i << "," << iterations;
				auto begin = std::chrono::high_resolution_clock::now();
				auto end = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double, std::milli> time;
				begin = std::chrono::high_resolution_clock::now();
				//render.RenderMandelMulti(xr1, xr2, yr1, yr2, iterations);







				MPI_Status stat;
				int world_size;
				MPI_Comm_size(MPI_COMM_WORLD, &world_size);
				ControllerEqualRange con(h, w, world_size - 1);
				//con.print();
				int** pixels = new int* [w];
				for (int i = 0; i < w; i++) {
					pixels[i] = new int[h];
				}
				double _yStep = 1.0 * (abs(yr1 - yr2) / (double)h);
				double _xStep = 1.0 * (abs(xr1 - xr2) / (double)w);
				//cout << "NEW iteration ###############################################" << std::endl;
				while (!con.Finished()) {
					if (sourceWaiting.size() == world_size - 1) {
						//cout << "cores ready: ";
						//for (auto& it : sourceWaiting) {
						//	cout << it << ", ";
						//}
						//cout << std::endl;
						for (auto& coreID : sourceWaiting) {
							//cout << "Resend new work to slave: " << coreID << std::endl;
							RangeParameters* param = con.GetWork(coreID);
							double msg_buffer[] = { xr1,yr1,static_cast<double>(param->x),static_cast<double>(param->y),static_cast<double>(h),_xStep, _yStep, static_cast<double>(maxIters) };
							MPI_Send(msg_buffer, 8, MPI_DOUBLE, coreID, TAG_INIT, MPI_COMM_WORLD);
						}
						sourceWaiting.clear();
					}
					else {
						MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
						int slave_rank = stat.MPI_SOURCE;
						// Decide according to the tag which type of message we have got
						if (stat.MPI_TAG == TAG_INIT) {
							//cout << "MASTER INIT from slave:" << slave_rank << std::endl;
							RangeParameters* param = con.GetWork(slave_rank);
							//no more jobs
							double* buff = new double[0];
							MPI_Recv(buff, 0, MPI_DOUBLE, slave_rank, TAG_INIT, MPI_COMM_WORLD, &stat);
							if (param == NULL) {
								//cout << "MASTER NO NEW JOB FOR SLAVE: " << slave_rank << std::endl;
								MPI_Send(NULL, 0, MPI_DOUBLE, slave_rank, TAG_NO_NEW_JOB, MPI_COMM_WORLD);
							}
							//Sends true init
							else {
								double msg_buffer[] = { xr1,yr1,static_cast<double>(param->x),static_cast<double>(param->y),static_cast<double>(h),_xStep, _yStep, static_cast<double>(maxIters) };
								MPI_Send(msg_buffer, 8, MPI_DOUBLE, slave_rank, TAG_INIT, MPI_COMM_WORLD);
							}
						}
						else if (stat.MPI_TAG == TAG_UPDATE) {
							
							double* buff = new double[updateWsize * h + 2];
							MPI_Recv(buff, updateWsize * h + 2, MPI_DOUBLE, slave_rank, TAG_UPDATE, MPI_COMM_WORLD, &stat);
							int newW = static_cast<int> (buff[0]);
							int newWend = static_cast<int> (buff[1]);
							//cout << "MASTER UPDATE from slave:" << slave_rank<<"  Range: "<<newW<<"-"<<newWend << std::endl;
							delete[]buff;
							RangeParameters* param = con.GetWork(slave_rank);
							if (newWend < newW) {
								param->x = newWend;
							}
							else {
								param->x = newW;
							}
							
							double* msg_buffer = new double[2];
							msg_buffer[0] = static_cast<double>(param->x);
							msg_buffer[1] = static_cast<double>(param->y);
							MPI_Send(msg_buffer, 2, MPI_DOUBLE, slave_rank, TAG_UPDATE, MPI_COMM_WORLD);
							delete[]msg_buffer;
						}
						else if (stat.MPI_TAG == TAG_JOB_FINISHED) {
							//cout << "MASTER JOB_FINISHED from slave:" << slave_rank << std::endl;
							double* buff = new double[updateWsize * h + 2];
							MPI_Recv(buff, updateWsize * h + 2, MPI_DOUBLE, slave_rank, TAG_JOB_FINISHED, MPI_COMM_WORLD, &stat);
							//cout << "MASTER JOB_FINISEHD RECV from slave:" << slave_rank << std::endl;
							int newW = static_cast<int> (buff[0]);
							int newWend = static_cast<int> (buff[1]);
							delete[]buff;
							RangeParameters* param = con.GetWork(slave_rank);
							if (newW > newWend) {
								param->x = newWend;
							}
							else {
								param->x = newW;
							}
							
							//con.print();
							//cout << "MASTER GETTING NEW JOB FOR SLAVE: " << slave_rank << std::endl;
							if (con.Finished()) {
								double* msg_buffer = NULL;
								sourceWaiting.push_back(slave_rank);
								MPI_Send(NULL, 0, MPI_DOUBLE, slave_rank, TAG_NO_NEW_JOB, MPI_COMM_WORLD);
								//cout << "MASTER TAG_NO_NEW_JOB TO slave:" << slave_rank << std::endl;

								//cout << "MASTER RENDER FINISHED____________________________________________________" << std::endl;
								/*cout << "cores ready: ";
								for (auto& it : sourceWaiting) {
									cout << it << ", ";
								}
								cout << std::endl;
								con.print();*/
								break;
							}
							param = con.GetWork(slave_rank);
							//cout << "MASTER RANGE of slave: " << slave_rank << "  range:" << newW<<"-"<<newWend<<std::endl;
							if (param != NULL) {
								if (param->y - param->x > 0) {
									double msg_buffer[] = { xr1,yr1,static_cast<double>(param->x),static_cast<double>(param->y),static_cast<double>(h),_xStep, _yStep, static_cast<double>(maxIters) };
									MPI_Send(msg_buffer, 8, MPI_DOUBLE, slave_rank, TAG_NEW_JOB, MPI_COMM_WORLD);
									//cout << "MASTER TAG_NEW_JOB TO slave:" << slave_rank << std::endl;
									//cout << "\tRange: " << param->x << "-" << param->y <<" True range: "<< param->y - param->x << std::endl;
								}
								else {
									double* msg_buffer = NULL;
									sourceWaiting.push_back(slave_rank);
									MPI_Send(NULL, 0, MPI_DOUBLE, slave_rank, TAG_NO_NEW_JOB, MPI_COMM_WORLD);
									//cout << "MASTER TAG_NO_NEW_JOB TO slave:" << slave_rank << std::endl;
								}

							}
							else {
								//cout << "\tcores ready: ";
								//for (auto& it : sourceWaiting) {
								//	cout << it << ", ";
								//}
								//cout << std::endl;
								sourceWaiting.push_back(slave_rank);
								double* msg_buffer = new double{0};
								MPI_Send(NULL, 0, MPI_DOUBLE, slave_rank, TAG_NO_NEW_JOB, MPI_COMM_WORLD);
								delete msg_buffer;
								//cout << "MASTER TAG_NO_NEW_JOB TO slave:" << slave_rank << std::endl;
							}
						}
						else {
							//cout << "MASTER BAD_TAG from slave: " << slave_rank << std::endl;
							//cout << "\t Error: " << stat.MPI_ERROR << std::endl;
							//cout << "\t TAG: " << stat.MPI_TAG << std::endl;
						}
					}
					/*if (i == 10 && iterations == 200) {
						cout << "\tcores ready: ";
						for (auto& it : sourceWaiting) {
							cout << it << ", ";
						}
						cout << std::endl;
					}*/
				}
				//cout << "MASTER FINISHED" << std::endl;
				end = std::chrono::high_resolution_clock::now();
				time = end - begin;
				if (rank == 0)
					cout << "," << time.count() << std::endl;
				for (int i = 0; i < w; i++) {
					delete[]pixels[i];
				}
				delete[]pixels;
				
			}
		}
	}
};
void Master() {
	//test1 64*64
	
	int ww = 64;
	int wh = 64;
	Window w(wh, ww);
	cout << "RESOLUTION: 64x64" << std::endl;
	w.RunTest();

	/*ww = 256;
	wh = ww;
	w.changeRes(wh, ww);
	cout << "RESOLUTION: 256x256" << std::endl;

	w.RunTest();

	ww = 512;
	wh = ww;
	w.changeRes(wh, ww);
	cout << "RESOLUTION: 512x512" << std::endl;
	w.RunTest();*/

	/*ww = 1024;
	wh = ww;
	w.changeRes(wh, ww);
	cout << "RESOLUTION: 1024x1024" << std::endl;
	w.RunTest();*/

	ww = 2048;
	wh = ww;
	w.changeRes(wh, ww);
	cout << "RESOLUTION: 2048x2048" << std::endl;
	w.RunTest();

	w.close();

}

void Slave() {

	MPI_Status stat;
	double offsetX;
	double offsetY;
	int wStart;
	int wEnd;
	int h;
	double wStep;
	double yStep;
	double maxIters;
	double* pixels = NULL;
	double* msg_buffer;
	// Here we send a message to the master asking for a job
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	bool first = true;
	while (true) {
		if (pixels == NULL) {
			if (first) {
				first = false;
				MPI_Send(NULL, 0, MPI_DOUBLE, 0, TAG_INIT, MPI_COMM_WORLD);
				MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
				if (stat.MPI_TAG == TAG_INIT) {
					// Retrieve job data from master into msg_buffer
					// [offsetX, offsetY, w.Start, w.end, h, wStep, hStep, maxIters
					msg_buffer = new double[8];
					MPI_Recv(msg_buffer, 8, MPI_DOUBLE, 0, TAG_INIT, MPI_COMM_WORLD, &stat);
					offsetX = msg_buffer[0];
					offsetY = msg_buffer[1];
					wStart = static_cast<int>(msg_buffer[2]);
					wEnd = static_cast<int>(msg_buffer[3]);
					h = static_cast<int>(msg_buffer[4]);
					wStep = msg_buffer[5];
					yStep = msg_buffer[6];
					maxIters = static_cast<double>(msg_buffer[7]);
					pixels = new double[updateWsize * h + 2];
					delete[]msg_buffer;
					//cout << "Good 1 INIT -- RANK:" << rank << "  Range:" << wStart << "-" << wEnd << std::endl;
				}
			}
			else {
				MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
				//cout << "2 Good 1 Probe -- RANK:" << rank << std::endl;
				if (stat.MPI_TAG == TAG_INIT) {
					// Retrieve job data from master into msg_buffer
					// [offsetX, offsetY, w.Start, w.end, h, wStep, hStep, maxIters
					msg_buffer = new double[8];
					MPI_Recv(msg_buffer, 8, MPI_DOUBLE, 0, TAG_INIT, MPI_COMM_WORLD, &stat);
					offsetX = msg_buffer[0];
					offsetY = msg_buffer[1];
					wStart = static_cast<int>(msg_buffer[2]);
					wEnd = static_cast<int>(msg_buffer[3]);
					h = static_cast<int>(msg_buffer[4]);
					wStep = msg_buffer[5];
					yStep = msg_buffer[6];
					maxIters = static_cast<double>(msg_buffer[7]);
					pixels = new double[updateWsize * h + 2];
					delete[]msg_buffer;
					//cout << "2 Good 1 INIT -- RANK:" << rank << "  Range:" << wStart << "-" << wEnd << std::endl;
				}
				else if (stat.MPI_TAG == TAG_TERMINATE) {
					
					MPI_Recv(NULL, 0, MPI_DOUBLE, 0, TAG_TERMINATE, MPI_COMM_WORLD, &stat);
					//cout << "Good 1 TERMINATE -- RANK:" << rank << std::endl;
					pixels = NULL;
					//break;
				}
				else if (stat.MPI_TAG == TAG_NO_NEW_JOB) {
					double* buff = new double[0];
					MPI_Recv(NULL, 0, MPI_DOUBLE, 0, TAG_INIT, MPI_COMM_WORLD, &stat);
					delete[] buff;
				}
				else {
					//bad tag
				}

			}
		}
		
		if (pixels == NULL) {
			//cout << "REAL TERMINATE IN SLAVE:" << rank << std::endl;
			break;
		}
		//loop
		int sendCounter = 0;
		for (; wStart < wEnd; wStart++) {

			for (int i = 0; i < h; i++) {
				if (updateWsize * sendCounter + h + 2 >= updateWsize*h+2) {
					//cout << "SEGMENTATION FAULT HERE" << std::endl;
				}
				pixels[updateWsize * sendCounter + h + 2] = calcMandelFracValue((1.0 * ((double)wStart)) * wStep + offsetX, ((double)i) * yStep + offsetY, maxIters);
			}
			sendCounter += 1;
			if (sendCounter == updateWsize-1) {
				pixels[0] = static_cast<double>(wStart);
				pixels[1] = static_cast<double>(wEnd);
				//cout << "SEND UPDATE SLAVE: " << rank<< "  Range:" << wStart << "-" << wEnd << std::endl;
				MPI_Send(pixels, updateWsize * h + 2, MPI_DOUBLE, 0, TAG_UPDATE, MPI_COMM_WORLD);
				
				MPI_Probe(0, TAG_UPDATE, MPI_COMM_WORLD, &stat);
				msg_buffer = new double[2];
				if (stat.MPI_TAG == TAG_UPDATE) {
					MPI_Recv(msg_buffer, 2, MPI_DOUBLE, 0, TAG_UPDATE, MPI_COMM_WORLD, &stat);
					//cout << "Good 2 RECV -- RANK:" << rank << std::endl;
				}
				wStart = static_cast<int>(msg_buffer[0]);
				wEnd = static_cast<int>(msg_buffer[1]);
				delete[]msg_buffer;
				if (wStart >= wEnd) {
					break;
				}
				sendCounter = 0;
			}
		}
		if (wStart > wEnd) {
			pixels[0] = wEnd;
			pixels[1] = wEnd;
		}
		else {
			pixels[0] = wStart;
			pixels[1] = wEnd;
		}
		
		//cout << "JOB_FINISHED SLAVE: " << rank << "   Range:" << wStart << "-" << wEnd << std::endl;
		MPI_Send(pixels, updateWsize * h + 2, MPI_DOUBLE, 0, TAG_JOB_FINISHED, MPI_COMM_WORLD);
		//cout << "Good 3 SEND -- RANK:" << rank << std::endl;
		
		MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
		delete[]pixels;
		pixels = NULL;

		if (stat.MPI_TAG == TAG_NEW_JOB) {
			// [offsetX, offsetY, w.Start, w.end, h, wStep, hStep, maxIters
			msg_buffer = new double[8];
			MPI_Recv(msg_buffer, 8, MPI_DOUBLE, 0, TAG_NEW_JOB, MPI_COMM_WORLD, &stat);
			//cout << "Good 3 TAG_NEW_JOB -- RANK:" << rank << std::endl;
			offsetX = msg_buffer[0];
			offsetY = msg_buffer[1];
			wStart = static_cast<int>(msg_buffer[2]);
			wEnd = static_cast<int>(msg_buffer[3]);
			h = static_cast<int>(msg_buffer[4]);
			wStep = msg_buffer[5];
			yStep = msg_buffer[6];
			maxIters = static_cast<double>(msg_buffer[7]);
			pixels = new double[updateWsize * h + 2];
			delete[]msg_buffer;
		}
		if (stat.MPI_TAG == TAG_NO_NEW_JOB) {
			MPI_Recv(NULL, 0, MPI_DOUBLE, 0, TAG_NO_NEW_JOB, MPI_COMM_WORLD, &stat);
			//cout << "Good 3 TAG_NO_NEW_JOB -- RANK:" << rank << std::endl;
			pixels = NULL;
		}
	}
}
int main(int argc, char** argv)
{
	argc_ = argc;
	argv_ = argv;
	MPI_Init(&argc_, &argv_);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		Master();
	}
	else {
		Slave();
	}

	MPI_Finalize();
	return 0;
}

