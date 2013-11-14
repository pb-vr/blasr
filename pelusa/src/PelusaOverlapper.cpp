#include "PelusaOverlapper.h"

// we have to use global variables until we can switch
// over to a more elegant thread solution
/*
 * RefSeq          *g_ref;
ReadClass       *g_read_a;
ifstream        *g_fin_a;
pthread_mutex_t g_mutex_fin=PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t g_mutex_fout=PTHREAD_MUTEX_INITIALIZER;
// JMS
ref_id_t        g_longReadUid = 0;
bit32_t         *g_n_aligned;
*/

PelusaOverlapper::PelusaOverlapper() 
	: 	debug(false),
		queryFile(""),
		targetFile(""),
		numProcs(1),
		kmerLength(8),
		encoder(NULL)
{
	
}

PelusaOverlapper::~PelusaOverlapper()
{
	for (int featureIdx = 0; featureIdx < numFeatures; featureIdx++)
	{
		delete(blooms[featureIdx]);
	}
	delete(encoder);
	/*TODO  take out completely for(map<pair<uint,uint>, vector<string>* >::const_iterator it = hashToIds.begin(); it != hashToIds.end(); ++it)
	{
		delete( (*it).second );
	}*/
	delete(byte2shorts);
}

void PelusaOverlapper::run()
{
    /* TODO remove
     * below is benchmarking code
     * int startTime = time(0);
    cerr << "Allocating GB array " << time(0) - startTime << endl;
    int array_size = 1073741824;
    int *big_array = new(int[array_size]);
    cerr << "Populating GB array " << time(0) - startTime << endl;
    for (int idx=0; idx < array_size; idx++)
    {
        big_array[idx] = idx;
    }
    cerr << "Finished populating array " << time(0) - startTime <<  endl;
    int sub_array_size = 65536;
    int *sub_array = new(int[sub_array_size]);
    cerr << "Summing 16384 64k blocks. " << time(0) - startTime << endl;
    for (int subIdx=0; subIdx < 16384; subIdx++)
    {
        // (rand() % 16384) * sub_array_size;
        // int startIdx = subIdx * sub_array_size; 
        int startIdx = (rand() % 16384) * sub_array_size;
        for (int sumIdx =0 ; sumIdx < sub_array_size; sumIdx++)
        {
            sub_array[sumIdx] = big_array[startIdx + sumIdx];
        }
    }
    // changing to sum above adds a second
    cerr << "Finished summing 16384 64k blocks. " << 
            time(0) - startTime << endl;
            */

	encoder = new FeatureEncoder(kmerLength);
	numFeatures = encoder->getNumFeatures(); // TODO should be in constructor?
    cerr << "Initializing bloom filters" << endl; // TODO put in real logging
	initializeBlooms();
    cerr << "Populating bloom filters" << endl; // TODO put in real logging
	populateBlooms();
    cerr << "Querying bloom filters" << endl; // TODO put in real logging
	queryBlooms();
    cerr << "Finished pelusa" << endl; // TODO put in real logging
}

void PelusaOverlapper::initializeBlooms()
{
	// create our bloom filters
	for (int featureIdx = 0; featureIdx < numFeatures; featureIdx++)
	{
		blooms.push_back( new bit_array_c( bloomWidth ) );
	}
	
	// populate the byte2ushortptr array which "spreads out" bits in a byte
	// into an eight byte unsigned long long int 
	// this data structure is useful for rapidly summing blooms in the query stage
	int bitsInByte = 8;
	byte2shorts = new BYTE_TO_SHORT_TYPE(boost::extents[256][8]);
	for (int shortIdx = 0; shortIdx < PELUSA_BYTE_TO_SHORT; shortIdx++)
	{	
		for (int bitIdx = 0; bitIdx < bitsInByte; bitIdx++)
		{	
			// TODO make sure we understand 7 - below
			// 7 -  below because we want the zeros place to be in the least significant bit.
			(*byte2shorts)[7 - bitIdx][bitIdx] = (shortIdx >> bitIdx) & 1;
		}
		//cout << byte2ushortptr[byteIdx] << endl;
	}
	
}

void PelusaOverlapper::populateBlooms()
{     
	FILE* file = fopen(targetFile.c_str(), "r");
    int recordCount = 0;
	while(true)
	{
		FastaRecord * record = new FastaRecord();
        if (++recordCount % 100 == 0)
        {
            cerr << "Populated record " << recordCount << endl;
        }
		if (!record->parseRecord(file))
		{	
			delete(record);
			break;
		}
		if (record->sequence.length() < (uint) kmerLength)
		{
			cerr << "Sequence of length " << record->sequence.length() << " is too small for k of length " << kmerLength << endl;
			delete(record);
			continue;
		}
		addRecordFeatures(record);	
		
		delete(record);
	}
	fclose(file);
    
    // if (debug) cerr << this->toString() << endl; 
    cerr << "Pelusa collision count " << this->collisionCount << endl;
} 


// overload for different hash functions, number etc.
void PelusaOverlapper::addRecordFeatures(FastaRecord * record)
{	
	// calculate hash values
	uint firstHashIdx  = RSHash(record->name) % bloomWidth;
	// uint secondHashIdx = JSHash(record->name) % bloomWidth;
	uint secondHashIdx = JSHash(record->name) % bloomWidth;
    // TODO pull out into routine for use below as well
	pair<uint, uint> key = 	firstHashIdx < secondHashIdx 				?
		 				 	pair<uint, uint>(firstHashIdx, secondHashIdx) 	:
			 				pair<uint, uint>(secondHashIdx, firstHashIdx);
	
	if (debug) 
        cerr << "Hash id for " << record->name << " = " << key.first << " - " << key.second << endl; 
	
	// add to our global map for mapping hash values back to record names

    while(true)
    {
        if(hashToIds.count(key) == 0)
            break;
        collisionCount += 1;
        secondHashIdx += collisionCount;
        secondHashIdx %= bloomWidth;
        key = firstHashIdx < secondHashIdx 				?
		 	    pair<uint, uint>(firstHashIdx, secondHashIdx) 	:
			 	pair<uint, uint>(secondHashIdx, firstHashIdx);
    }
	hashToIds.insert( pair < pair<uint, uint>, string > (key, record->name));
	
	// set the bits at the hashed index positions in each feature's bloom filter
	// TODO replace with a stack array for efficiency?
	vector<uint>* features = new vector<uint>; 
	encoder->encode(record->sequence, features);
	for (uint featureIdx=0; featureIdx < features->size(); featureIdx++)
	{
		// if (debug) cerr << "Feature " << featureIdx << "=" << (*features)[featureIdx] << endl;
		bit_array_c* featureBloom = blooms[ (uint)(*features)[featureIdx] ];
		featureBloom->SetBit(firstHashIdx);
		featureBloom->SetBit(secondHashIdx);	
	}
	
	delete(features);
}

void* threadStarter(void* arg)
{
	void * dummy;
    // TODO better return value
	((PelusaWorker*)arg)->run();	
	return dummy;
}

void PelusaOverlapper::queryBlooms()
{     
	FILE* file = fopen(queryFile.c_str(), "r");
	pthread_mutex_t queryMutex=PTHREAD_MUTEX_INITIALIZER;
	pthread_mutex_t outputMutex=PTHREAD_MUTEX_INITIALIZER;
	
	vector<pthread_t> threads(numProcs);
	vector<PelusaWorker *> workers(numProcs);
	// create threads
	for (int workerIdx=0; workerIdx < numProcs; workerIdx++)
	{
		cerr << "Creating worker " << workerIdx << endl;
		PelusaWorker * worker = new PelusaWorker(workerIdx);
		workers[workerIdx] = worker;
		worker->setOverlapper(this);
		worker->setQueryFh(file);
		worker->setQueryMutexPtr(&queryMutex);
		worker->setOutputMutexPtr(&outputMutex);
		pthread_create(&threads[workerIdx], NULL, threadStarter, (void*)worker);
	}
	
	for (int workerIdx=0; workerIdx < numProcs; workerIdx++)
	{
		pthread_join(threads[workerIdx], NULL);
	}
	
	for (int workerIdx=0; workerIdx < numProcs; workerIdx++)
	{
		delete(workers[workerIdx]);
	}
	
	fclose(file);
}

int PelusaOverlapper::queryRecord(
        FastaRecord * record, 
        map<string, int>* id2score)
{
	if (record->sequence.length() < (uint) kmerLength)
	{
		cerr << "Sequence of length " << record->sequence.length() << 
             " is too small for k of length " << kmerLength << endl;
		return 0;
	}
	vector<uint>* features = new vector<uint>; 
	encoder->encode(record->sequence, features);
	
  	int* sumFeatures = new(int[bloomWidth]);
    std::fill(sumFeatures, sumFeatures + bloomWidth, 0);
    
    for (uint featureIdx = 0; featureIdx < features->size(); featureIdx++)
    {
    	uint bloomIdx = (uint)(*features)[featureIdx];
        for (int bitIdx = 0; bitIdx < bloomWidth; bitIdx++)
        {  			      
            sumFeatures[bitIdx] += (*blooms[bloomIdx])[bitIdx];
        }
    }

    if (debug)
    {
        cerr << "Sum vector for query " << record->name << endl;
        for (int bitIdx = 0; bitIdx < bloomWidth; bitIdx++)
        {  			      
            cerr << bitIdx << ": " << sumFeatures[bitIdx] << ", ";
        }
        cerr << endl;
    }
    
	// determine the maximum indices summed features array 
    vector<int>* topIndices = new(vector<int>);
    findTopIndices(sumFeatures, topIndices);
    
    // iterate through pairs of the valueIndex and pull out pairs of indices
    // that actually map to real ids. give the sum of the value as their scores.
    //  		for ( firstPairIterator=maximizer->begin() ; firstPairIterator != maximizer->end(); firstPairIterator++ )

    int numHits = 0;
    for (int idx=0; idx < topIndices->size(); idx++)
    {
        int topIndexI = (*topIndices)[idx];
        for (int jdx=idx+1; jdx < topIndices->size(); jdx++)
        {
            int topIndexJ = (*topIndices)[jdx];
        	pair<uint, uint> key = 	topIndexI < topIndexJ ?
		 				 	pair<uint, uint>(topIndexI, topIndexJ) 	:
			 				pair<uint, uint>(topIndexJ, topIndexI);

            int score = min(sumFeatures[topIndexI], sumFeatures[topIndexJ]);
            if (debug)
            {
                cerr << "Score for " << topIndexI << "," << 
                                        topIndexJ << ":" << 
                                        score << endl;
            }
            string hit;
            if (hashToIds.count(key) != 0) 
            {
                hit = hashToIds.find(key)->second;
                id2score->insert(pair<string, int>(hit, score));
                numHits++;
            }
        }
    }
	
	return numHits;
}

// populates topIndices vector with the "topColumns" top indices from sumFeatures
void PelusaOverlapper::findTopIndices(int* array, vector<int>* topIndices )
{
    list< pair<int,int> >* valueIndexPairs = new list< pair<int,int> >;
    for (int idx = 0; idx < bloomWidth; idx++)
    {
        pair<int, int> newPair = pair<int,int>(array[idx],idx);
        valueIndexPairs->push_back(newPair);
    }
    valueIndexPairs->sort();
    valueIndexPairs->reverse();
    valueIndexPairs->resize(topColumns);
    list< pair<int, int> >::iterator pairIterator;
    for (pairIterator =  valueIndexPairs->begin(); 
         pairIterator != valueIndexPairs->end(); 
         pairIterator++)
    {
        topIndices->push_back((*pairIterator).second);
    }
    return;
}


void PelusaOverlapper::setBloomWidth(int bloomWidth)
{
	if ( floor(bloomWidth / 8) != bloomWidth / 8 )
    {
    	throw PelusaException("Bloom width should be divisible by 8!");
    }
	this->bloomWidth = bloomWidth;
}


const string PelusaOverlapper::toString()
{
	stringstream stream;
	stream << "queryFile   = " << queryFile   << endl;
	stream << "targetFile  = " << targetFile  << endl;
	stream << "numProcs    = " << numProcs    << endl;
	stream << "kmerLength  = " << kmerLength  << endl;		
	stream << "bloomWidth  = " << bloomWidth  << endl;	
	stream << "numSegments = " << numSegments << endl;
	stream << "topColumns  = " << topColumns  << endl;
	for (int bloomIdx = 0; bloomIdx < numFeatures; bloomIdx++)
	{	
		stream << bloomIdx << ": ";
		blooms[bloomIdx]->Dump(stream);
		stream << endl;
	}
	
	return stream.str();
	
}


PelusaWorker::PelusaWorker(int id)
	:	id(id)
{
	
}


void PelusaWorker::run()
{
	cerr << "Starting worker " << id << endl;
    int recordCount = 0;
	while(true)
	{
		FastaRecord * record = new FastaRecord();
		
		// get the input
		pthread_mutex_lock(queryMutexPtr);
		bool wasParsed = record->parseRecord(queryFh)	;
        recordCount++;
		pthread_mutex_unlock(queryMutexPtr);
		
		// termination condition for all threads
		if (!wasParsed)
		{	
			delete(record);
			break;
		}
		
        if (recordCount % 100 == 0)
        {
            cerr << "Querying record " << recordCount 
                 << " in worker " << id << endl;
        }

		// the heavy lifting
		map<string, int>* id2score = new map<string, int>;
		// int numHits = overlapper->queryRecordPaired(record, id2score);
		int numHits = overlapper->queryRecord(record, id2score);

		// output
		if (numHits > 0){
			pthread_mutex_lock(outputMutexPtr);
			outputScores(record, id2score);
			pthread_mutex_unlock(outputMutexPtr);
		}
		
		// cleanup
		delete(record);
		delete(id2score); // TODO put both of these on stack?
	}
	cerr << "Finishing worker " << id << endl;
}

void PelusaWorker::setOverlapper(PelusaOverlapper * overlapper)
{
	this->overlapper = overlapper;
}

void PelusaWorker::setQueryFh(FILE* queryFh)
{
	this->queryFh = queryFh;
}

void PelusaWorker::setQueryMutexPtr(pthread_mutex_t * queryMutexPtr)
{
	this->queryMutexPtr = queryMutexPtr;	
}

void PelusaWorker::setOutputMutexPtr(pthread_mutex_t * outputMutexPtr)
{
	this->outputMutexPtr = outputMutexPtr;	
}

// TODO check query parsing parellization. I think it's okay serially, but?

// TODO put in filehandle argument?
void PelusaWorker::outputScores(FastaRecord * record, map<string, int>* id2score)
{
	for(map<string, int>::const_iterator it = id2score->begin(); it != id2score->end(); ++it)
	{
		cout << record->name << "\t" << (*it).first << "\t" << (*it).second << "\n"; 
	}
}
