#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include <cstdint>
#include <boost/filesystem.hpp>

#include "CImg.h"

using namespace std;
using namespace boost::filesystem;
using namespace cimg_library;



vector<string> GetImageList(string dir)
{
	path p(dir);

	directory_iterator end_itr;

	vector<string> retFiles;
	static const vector<string> imgExtensions = { ".jpg",".png",".bmp" };

	// cycle through the directory
	for (auto itr: recursive_directory_iterator(p))
	{
		// If it's not a directory, list it. If you want to list directories too, just remove this check.
		
		if (is_regular_file(itr.path())) {
			path ext = itr.path().extension();
			if (find(imgExtensions.begin(),imgExtensions.end(),ext.string())!=imgExtensions.end())
				retFiles.push_back(itr.path().string());
		}
	}
	return retFiles;
}

class ImgBlock
{
public:
	ImgBlock() :clr_{ 0,0,0 } {}
	ImgBlock(CImg<unsigned char> img, string file) :img_(img), file_(file){}
	ImgBlock(const ImgBlock& rhs)
	{
		img_ = rhs.img_;
		clr_[0] = rhs.clr_[0];
		clr_[1] = rhs.clr_[1];
		clr_[2] = rhs.clr_[2];
		file_ = rhs.file_;
	}

	void setAvgColor(unsigned char r, unsigned char g, unsigned char b)
	{
		clr_[0] = r;
		clr_[1] = g;
		clr_[2] = b;
	}

	const unsigned char* color() const
	{
		return clr_;
	}

	string file() const
	{
		return file_;
	}

	CImg<unsigned char>& image()
	{
		return img_;
	}

	bool operator<(const ImgBlock& rhs) const
	{
		return (clr_[0]*clr_[0] + clr_[1]*clr_[1] + clr_[2]*clr_[2]) < (rhs.clr_[0]*rhs.clr_[0] + rhs.clr_[1]*rhs.clr_[1] + rhs.clr_[2]*rhs.clr_[2]);
	}

	//for searching, want to match file name
	bool operator==(const string &rhsFile) const
	{
		return file_ == rhsFile;
	}

private:
	CImg<unsigned char> img_;
	unsigned char clr_[3];
	string file_;
};


double cDiff(const unsigned char* clr1, const unsigned char* clr2)
{
	unsigned rAvg = (clr1[1] + clr2[1]) / 2;
	unsigned dR = clr1[1] - clr2[1];
	unsigned dG = clr1[2] - clr2[2];
	unsigned dB = clr1[3] - clr2[3];
	return(sqrt((2 + rAvg / 256)*dR*dR + 4 * dG*dG + (2 + (255 - rAvg) / 256)*dB*dB));
}

string findClosestShade(vector<ImgBlock>& blocks, unsigned char* clr, unsigned char* retAvg)
{
	double minDiff = cDiff(blocks.begin()->color(), clr);
	ImgBlock minDiffElement(*(blocks.begin()));

	for (vector<ImgBlock>::const_iterator it(blocks.begin()); it != blocks.end(); ++it)
	{
		double diff = cDiff(it->color(), clr);
		if (diff < minDiff)
		{
			minDiff = diff;
			minDiffElement = *it;
		}
	}
	retAvg[0] = minDiffElement.color()[0];
	retAvg[1] = minDiffElement.color()[1];
	retAvg[2] = minDiffElement.color()[2];
	return minDiffElement.file();
}

void applyTint(CImg<unsigned char>& img, double f)
{
	for (int x = 0; x < img.width(); x++)
	{
		for (int y = 0; y < img.height(); y++)
		{
			img(x, y, 0, 0) = img(x, y, 0, 0) + (255 - img(x, y, 0, 0)) * f;
			img(x, y, 0, 1) = img(x, y, 0, 1) + (255 - img(x, y, 0, 1)) * f;
			img(x, y, 0, 2) = img(x, y, 0, 2) + (255 - img(x, y, 0, 2)) * f;
		}
	}
}

void applyShade(CImg<unsigned char>& img, double f)
{
	for (int x = 0; x < img.width(); x++)
	{
		for (int y = 0; y < img.height(); y++)
		{
			img(x, y, 0, 0) = img(x, y, 0, 0) * (1 - f);
			img(x, y, 0, 1) = img(x, y, 0, 1) * (1 - f);
			img(x, y, 0, 2) = img(x, y, 0, 2) * (1 - f);
		}
	}
}

int main() {

	vector<string> imgList = GetImageList("D:\\Pictures");
	vector<ImgBlock> blocks;
	
	string inputImgFile;

	cout << "Enter mosaic image: ";
	cin >> inputImgFile;

	int outputWidth;
	cout << "Output image width: ";
	cin >> outputWidth;

	int blockSize;
	cout << "Block size: ";
	cin >> blockSize;

	for (vector<string>::const_iterator it(imgList.begin()); it != imgList.end(); ++it)
	{
		try
		{
			CImg<unsigned char> image(it->c_str());
			image.resize(blockSize, blockSize);
			ImgBlock newBlock(image, *it);
			image.resize(1, 1);
			unsigned char* pix = image.data();
			newBlock.setAvgColor(pix[0], pix[1], pix[2]);
			//pixel data is stored NOT interleaved, so for a normal image R1R2R3...G1G2G3...B1G2G3...
			//in this case it doesn't matter since I've resized it to be 1x1

			blocks.push_back(newBlock);
		}
		catch (...)
		{
			cout << "Error";
		}
			
	}



	CImg<unsigned char> inputImg(inputImgFile.c_str());
	int width = inputImg.width();
	int height = inputImg.height();
	inputImg.resize(outputWidth,outputWidth*height/width);
	width = inputImg.width();
	height = inputImg.height();
/*
	while (inputImg.width() % blockSize | inputImg.height() % blockSize)
	{
		blockSize++;
		if ((blockSize > inputImg.width()) | (blockSize > inputImg.height()))
		{
			cout << "Error: Try a smaller block size";
			return 0;
		}
	}*/

	
	//unsigned char* inputData = inputImg.data();


	vector<CImg<unsigned char> > outBlocks;

	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			cout << "Processing block (" << x << "," << y << ")" << endl;
			unsigned char inputPixelValue[3];
			inputPixelValue[0] = inputImg(x,y,0,0);
			inputPixelValue[1] = inputImg(x,y,0,1);
			inputPixelValue[2] = inputImg(x,y,0,2);
			//cout << "R: " << int(inputImg(x, y, 0, 0)) << " G: " << int(inputImg(x, y, 0, 1)) << " B: " << int(inputImg(x, y, 0, 2))<<endl;
			
			//find image with closest average pixel value
			unsigned char retAvg[3];
			string curBlockFile = findClosestShade(blocks, inputPixelValue, retAvg);

			
			//scale to block size
			vector<ImgBlock>::iterator curBlockIt = find(blocks.begin(), blocks.end(), curBlockFile);
			//blockImg.resize(blockSize, blockSize); //TODO: scale then crop to keep aspect intact
 		    /*
			//adjust brightness of image to match pixel
			int64_t inputDistSq = inputPixelValue[0] * inputPixelValue[0] + inputPixelValue[1] * inputPixelValue[1] + inputPixelValue[2] * inputPixelValue[2];
			int64_t outputDistSq = retAvg[0] * retAvg[0] + retAvg[1] * retAvg[1] + retAvg[2] * retAvg[2];
			int64_t c = outputDistSq - inputDistSq;
			if (inputDistSq > outputDistSq)
			{
				int64_t ar = (255 - retAvg[0]);
				int64_t ag = (255 - retAvg[1]);
				int64_t ab = (255 - retAvg[2]);
				int64_t a = ar*ar + ag*ag + ab*ab;
				int64_t b = 2 * (ar + ab + ag);
				cout << sqrt(b*b - 4*a*c) << endl;
				double f = (-sqrt(b*b - 4*a*c) - b) / (2 * a);
				applyTint(curBlockIt->image(), f/2);
			}
			else
			{
				int64_t a = outputDistSq;
				int64_t b = -2 * outputDistSq;
				double f = (-sqrt((b*b) - (4*a*c)) - b) / (2 * a);
				applyShade(curBlockIt->image(), f/2);
			}
			*/
			//save block
			outBlocks.push_back(curBlockIt->image());
			//cout << "R: "<<inputPixelValue[0]<< " G: " << inputPixelValue[1]<< " B: " << inputPixelValue[2]<<endl;
		}
	}

	int outputHeight = outputWidth*height / width;
	CImg<unsigned char> outputImage(outputWidth*blockSize, outputHeight*blockSize,1,3); //create output image with correct dimensions

	//iterator over blocks and build up output mosaic
	{ int i = 0;
	for (vector<CImg<unsigned char>>::iterator it(outBlocks.begin()); it != outBlocks.end(); ++it, ++i)
	{
		
		for (int x = (i%outputWidth)*blockSize, xi = 0; xi < blockSize; ++xi, ++x)
		{
			for (int y = int(i / outputWidth)*blockSize, yi = 0; yi < blockSize; ++yi, ++y)
			{
				outputImage(x, y, 0, 0) = (*it)(xi,yi,0,0); //r
				outputImage(x, y, 0, 1) = (*it)(xi,yi, 0, 1); //g
				outputImage(x, y, 0, 2) = (*it)(xi,yi, 0, 2); //b	 
			}
			cout << endl;
		}
	}
	}
	CImgDisplay main_disp(outputImage, "Test Mosaic");
	while (!main_disp.is_closed())
		main_disp.wait();

	/*
	CImg<unsigned char> image("images/milla.bmp"), visu(500, 400, 1, 3, 0);
	const unsigned char red[] = { 255,0,0 }, green[] = { 0,255,0 }, blue[] = { 0,0,255 };
	image.blur(2.5);
	CImgDisplay main_disp(image, "Click a point"), draw_disp(visu, "Intensity profile");
	while (!main_disp.is_closed() && !draw_disp.is_closed()) {
		main_disp.wait();
		if (main_disp.button() && main_disp.mouse_y() >= 0) {
			const int y = main_disp.mouse_y();
			visu.fill(0).draw_graph(image.get_crop(0, y, 0, 0, image.width() - 1, y, 0, 0), red, 1, 1, 0, 255, 0);
			visu.draw_graph(image.get_crop(0, y, 0, 1, image.width() - 1, y, 0, 1), green, 1, 1, 0, 255, 0);
			visu.draw_graph(image.get_crop(0, y, 0, 2, image.width() - 1, y, 0, 2), blue, 1, 1, 0, 255, 0).display(draw_disp);
		}
	}
	*/
	return 0;
}