#ifndef SDP_FRAGMENT_H_
#define SDP_FRAGMENT_H_


class Fragment {
 public:
	unsigned int x;
	unsigned int y;
	unsigned int weight;
	unsigned int length;
	int index;
	int chainPrev;
	int cost;
    int above;
	unsigned int chainLength;

	unsigned int GetX() const {
		return x;
	}
	unsigned int GetY() const {
		return y;
	}

	Fragment(unsigned int px, unsigned int py, int pweight=0) {
		x = px;
		y = py;
		weight = pweight;
        length = index = 0;
		chainPrev = cost = chainLength = 0;
        cost = 0;
        above = -1;
	}

    int SetAbove(int a) {
        above = a;
        return true;
    }

    bool GetAbove(int & a) {
        if (above >= 0) {
            a = above;
            return true;
        } else {
            a = -1;
            return false;
        }
    }

    //
    // Provide default constructor that will
    // give bad results if members are not properly initialized
    // later on.
    //
	Fragment() {
		x = -1;
		y = -1;
        weight = length = index = 0;
		chainPrev = cost = chainLength = 0;
        above = -1;
	}

	int LessThanXY(const Fragment &f) const {
        if (x == f.x) {
            if (y == f.y) {
                if (length == f.length) return 0;
                else return length < f.length;
            } else return y < f.y;
        } else return x < f.x;

        /*
		if (x < f.x)
			return 1;
		else if (x == f.x) 
			return y < f.y;
		else 
			return 0;
        */
	}

	int LessThanYX(const Fragment &f) const {
        if (y == f.y) {
            if (x == f.x) {
                if (length == f.length) return 0;
                else return length < f.length;
            } else return x < f.x;
        } else return y < f.y;

        /*
		if (y < f.y)
			return 1;
		else if (y == f.y) 
			return x < f.x;
		else 
			return 0;
        */
	}

	int operator<(const Fragment &f) const {
		// 
		// Sort fragments by diagonal:
		//
		int diag, fDiag;
		diag = (y - x);
		fDiag = f.y - f.x;
		if (diag < fDiag)
			return 1;
		else if (diag == fDiag)
			return (x < f.x);
		else
			return 0;
	}
	Fragment& operator=(const Fragment &rhs) {
		x           = rhs.x;
		y           = rhs.y;
		index       = rhs.index;
		cost        = rhs.cost;
		weight      = rhs.weight;
        length      = rhs.length;
		chainLength = rhs.chainLength;
		chainPrev   = rhs.chainPrev;
        above       = rhs.above;
		return *this;
	}
		
	int operator==(const Fragment &f) const {
		return (x == f.x and y == f.y);
	}
	int operator>(const Fragment &f) const {
		return (!(*this < f) &&  !(*this == f));
	}
	int GetLength() {
		return length;
	}
	void SetLength(int _length) {
		length = _length;
	}
};


#endif
