#ifndef __TIMER__
#define __TIMER__

class Timer {
    private:
        // dim holds total days of month
        // also note that m_day is the cumulative day number in the simulation (not the day of the month!)
        const std::vector<int> dim { 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365}; // last yday in each month
        int m_step { 1 };
        int m_halfhour { 0 };
        int m_day { 1 };
        int m_yday { 1 }; 
        int m_month { 1 };
        int m_year { 1 };
        bool m_newday { true };
        bool m_newmonth { false }; 
        bool m_newyear { false };
        public:
            Timer(const int yday = 1) : m_yday(yday) {
                m_month = std::lower_bound(dim.begin(), dim.end(), yday) - dim.begin() + 1;
                if (m_month == 1) {
                    m_day = m_yday;
                } else {
                    m_day = m_yday - dim[m_month-2]; // subtract days of previous month(s)
                }
                if (m_day == 1) m_newmonth = true;
                if (m_yday == 1) m_newyear = true;
            };
        int next() { 
            m_newday = false;
            m_newmonth = false;
            m_newday = false;
            ++m_step;
            ++m_halfhour;
            if (m_halfhour > 48) {
                m_halfhour = 1;
                ++m_yday;
                ++m_day;
                m_newday = true;
                if (m_yday > dim[m_month-1]) {
                    if (++m_month > 12) {
                        ++m_year;
                        m_yday = 1;
                        m_month = 1;
                        m_newyear = true;
                    }
                    m_newmonth = true;
                }
            }
            return m_step;
        };
        const int step() { return m_step; }
        const int day() { return m_day; }
        const int yday() { return m_yday; }
        const int month() { return m_month; }
        const int quarter() { return ceil(m_month/3.0f); }
        const int year() { return m_year; }
        bool isNewYear() { return m_newyear; }
        bool isNewMonth() { return m_newmonth; }
        bool isNewDay() { return m_newday; }
        bool isNewQuarter() { return m_newmonth && ((m_month % 3) == 1); }
};

#endif // __TIMER__