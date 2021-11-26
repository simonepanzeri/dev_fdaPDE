//
// Created by simonepanzeri on 24/11/2021.
//

#ifndef DEV_FDAPDE_DE_DATA_IMP_H
#define DEV_FDAPDE_DE_DATA_IMP_H

template<UInt ndim>
DEData<ndim>::DEData(const std::vector<Point<ndim>>& data, const UInt& order, const VectorXr& fvec, Real heatStep,
                     UInt heatIter, const std::vector<Real>& lambda, const UInt& nfolds, const UInt& nsim,
                     const std::vector<Real>& stepProposals, Real tol1, Real tol2, bool print, UInt search) :
        data_(data), order_(order), fvec_(fvec), heatStep_(heatStep), heatIter_(heatIter), lambda_(lambda),
        Nfolds_(nfolds), nsim_(nsim), stepProposals_(stepProposals), tol1_(tol1), tol2_(tol2), print_(print),
        search_(search) {
}

template<UInt ndim>
void DEData<ndim>::printData(std::ostream & out) const
{
    for(int i = 0; i < data_.size(); i++)
    {
        out << data_[i] << std::endl;
    }
}

template <UInt ndim>
DEData_time<ndim>::DEData_time(const std::vector<Point<ndim>>& data, const std::vector<Real>& data_time) :
    data_time_(data_time) {
    static constexpr Real eps = std::numeric_limits<Real>::epsilon(), tolerance = 100 * eps;

    //extraction of non duplicated times
    std::set<Real> set_time_noD(data_time_.cbegin(),data_time_.cend());
    data_time_noD_.resize(set_time_noD.size());
    std::copy(set_time_noD.begin(), set_time_noD.end(), data_time_noD_.begin());
    //std::vector<Real> v(set_time_noD.begin(), set_time_noD.end());
    //data_time_noD_.swap(v);
    set_time_noD.clear();

    //creation of the map which contains as keys the IDs of different spatial points and as values the indices of
    //the time instants point in data_time_noD_ (this structure is needed to build one of the FEmatrices)
    std::set<UInt> set_helper;
    for(size_t i = data.front().id(); i <= data.back().id(); ++i) {
        if(!isAlready(i, set_helper)) {insert(i); set_helper.insert(i);}
        for(size_t j = i + 1; j <= data.back().id(); ++j){
            if (data[i].dist(data[j]) < tolerance and !isAlready(j, set_helper))
                {insert(i,j); set_helper.insert(j);}
        }
    }
    set_helper.clear();
}

template <UInt ndim>
const std::vector<UInt> DEData_time<ndim>::getID_noD() const {
    std::vector<UInt> noDup;
    noDup.reserve(data_noD_.size());
    for(const auto &s : data_noD_)
        noDup.push_back(s.first);
    return noDup;
}

template <UInt ndim>
const std::set<UInt>& DEData_time<ndim>::getTimesIndices(UInt id) const {
    return data_noD_.find(id)->second;
};

template <UInt ndim>
void DEData_time<ndim>::printData(std::ostream& out) const {
    for(int i = 0; i < data_time_.size(); i++)
    {
        out << data_time_[i] << std::endl;
    }
}

template<UInt ndim>
void DEData_time<ndim>::printMap(std::ostream& out) const {
    for(std::map<UInt, std::set<UInt>>::const_iterator it = data_noD_.cbegin(); it != data_noD_.cend(); ++it) {
        out << "id: " << it->first << '\t' << "times: ";
        for(const UInt& t : it->second)
            out << std::setw(3) << t;
        out << std::endl;
    }
}

#endif //DEV_FDAPDE_DE_DATA_IMP_H
