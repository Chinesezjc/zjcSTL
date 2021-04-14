// This Code was made by Chinese_zjc_.

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

namespace zjcSTL
{
    template <typename __Tp1, typename __Tp2>
    __Tp1 power(__Tp1 A, __Tp2 B)
    {
        __Tp1 res = 1;
        while (B > 0)
        {
            if (B % 2 == 1)
            {
                res *= A;
            }
            A *= A;
            B /= 2;
        }
        return res;
    }
    const double eps = 1e-10;
    const double INF = 1e50;
    template <class __Tp>
    class complex
    {
    public:
        __Tp x, y;
        _GLIBCXX14_CONSTEXPR complex(const __Tp &A = __Tp(), const __Tp &B = __Tp()) : x(A), y(B) {}
        complex operator+(const complex &__Val) const _GLIBCXX_NOEXCEPT { return complex(x + __Val.x, y + __Val.y); }
        complex operator-(const complex &__Val) const _GLIBCXX_NOEXCEPT { return complex(x - __Val.x, y - __Val.y); }
        complex operator*(const complex &__Val) const _GLIBCXX_NOEXCEPT
        {
            return complex(x * __Val.x - y * __Val.y, x * __Val.y + y * __Val.x);
        }
        complex operator/(const complex &__Val) const _GLIBCXX_NOEXCEPT
        {
            return complex((x * __Val.x + y * __Val.y) / (__Val.x * __Val.x - __Val.y * __Val.y),
                           (__Val.x * y - x * __Val.y) / (__Val.x * __Val.x - __Val.y * __Val.y));
        }
        complex &operator+=(const complex &__Val) _GLIBCXX_NOEXCEPT { return *this = *this + __Val; }
        complex &operator-=(const complex &__Val) _GLIBCXX_NOEXCEPT { return *this = *this - __Val; }
        complex &operator*=(const complex &__Val) _GLIBCXX_NOEXCEPT { return *this = *this * __Val; }
        complex &operator/=(const complex &__Val) _GLIBCXX_NOEXCEPT { return *this = *this / __Val; }
        friend complex conj(const complex &__Val) _GLIBCXX_NOEXCEPT { return complex(__Val.x, -__Val.y); }
        friend std::istream &operator>>(std::istream &ist, complex &__Val) { return ist >> __Val.x >> __Val.y; }
        friend std::ostream &operator<<(std::ostream &ost, const complex &__Val)
        {
            return ost << '(' << __Val.x << ',' << __Val.y << ')';
        }
    };
    template <typename __Tp>
    complex<__Tp> polar(const __Tp &__rho, const __Tp &__theta) _GLIBCXX_NOEXCEPT
    {
        __glibcxx_requires_subscript(__rho >= 0);
        return complex<__Tp>(__rho * cos(__theta), __rho * sin(__theta));
    }
    template <size_t base = 10>
    class bigfloat
    {
    private:
        std::vector<long long> number;
        bool mark;
        int point;
        int up_per_time;
        void shrink_to_fit() _GLIBCXX_NOEXCEPT
        {
            while (!(number.empty() || number.back()))
            {
                number.pop_back();
            }
            if (number.empty())
            {
                mark = false;
                point = 0;
            }
        }

    public:
        class iterator
        {
        private:
            std::vector<long long>::iterator Iter;

        public:
            _GLIBCXX14_CONSTEXPR iterator() _GLIBCXX_NOEXCEPT
            {
            }
            iterator(const std::vector<long long>::iterator &__Val) _GLIBCXX_NOEXCEPT
            {
                Iter = __Val;
            }
            iterator operator++() _GLIBCXX_NOEXCEPT
            {
                return ++Iter;
            }
            iterator operator++(int) _GLIBCXX_NOEXCEPT
            {
                return Iter++;
            }
            iterator operator--() _GLIBCXX_NOEXCEPT
            {
                return --Iter;
            }
            iterator operator--(int) _GLIBCXX_NOEXCEPT
            {
                return Iter--;
            }
            bool operator==(const iterator &__Val) const _GLIBCXX_NOEXCEPT
            {
                return Iter == __Val.Iter;
            }
            bool operator!=(const iterator &__Val) const _GLIBCXX_NOEXCEPT
            {
                return Iter != __Val.Iter;
            }
            long long &operator*() const _GLIBCXX_NOEXCEPT
            {
                return *Iter;
            }
        };
        class reverse_iterator
        {
        private:
            std::vector<long long>::reverse_iterator Iter;

        public:
            _GLIBCXX14_CONSTEXPR reverse_iterator() _GLIBCXX_NOEXCEPT
            {
            }
            reverse_iterator(const std::vector<long long>::reverse_iterator &__Val) _GLIBCXX_NOEXCEPT
            {
                Iter = __Val;
            }
            reverse_iterator operator++() _GLIBCXX_NOEXCEPT
            {
                return ++Iter;
            }
            reverse_iterator operator++(int) _GLIBCXX_NOEXCEPT
            {
                return Iter++;
            }
            reverse_iterator operator--() _GLIBCXX_NOEXCEPT
            {
                return --Iter;
            }
            reverse_iterator operator--(int) _GLIBCXX_NOEXCEPT
            {
                return Iter--;
            }
            bool operator==(const reverse_iterator &__Val) const _GLIBCXX_NOEXCEPT
            {
                return Iter == __Val.Iter;
            }
            bool operator!=(const reverse_iterator &__Val) const _GLIBCXX_NOEXCEPT
            {
                return Iter != __Val.Iter;
            }
            long long &operator*() const _GLIBCXX_NOEXCEPT
            {
                return *Iter;
            }
        };
        iterator begin() _GLIBCXX_NOEXCEPT
        {
            return number.begin();
        }
        iterator end() _GLIBCXX_NOEXCEPT
        {
            return number.end();
        }
        reverse_iterator rbegin() _GLIBCXX_NOEXCEPT
        {
            return number.rbegin();
        }
        reverse_iterator rend() _GLIBCXX_NOEXCEPT
        {
            return number.rend();
        }
        _GLIBCXX14_CONSTEXPR bigfloat(const int &Up_per_time = 10) _GLIBCXX_NOEXCEPT
            : number(),
              mark(),
              point(),
              up_per_time(Up_per_time) {}
        template <typename __Tp>
        bigfloat(const __Tp &__Val, const int Up_per_time = 10) _GLIBCXX_NOEXCEPT
            : number(),
              mark(),
              point(),
              up_per_time(Up_per_time)
        {
            __Tp tmp = __Val;
            if (tmp < 0)
            {
                mark = true;
                tmp = -tmp;
            }
            while (tmp - floor(tmp) > eps)
            {
                tmp *= base;
                --point;
            }
            // while()
        }
        bigfloat(const std::string &__Val, const int &Up_per_time = 10) _GLIBCXX_NOEXCEPT
        {
            number.clear();
            mark = false;
            point = 0;
            up_per_time = Up_per_time;
            std::string tmp = __Val;
            std::reverse(tmp.begin(), tmp.end());
            if (*--tmp.end() == '-')
            {
                mark = true;
                tmp.erase(--tmp.end());
            }
            if (base <= 36)
            {
                int i = 0;
                while (i < (int)tmp.length() && tmp[i] != '.')
                {
                    number.push_back(tmp[i] <= '9' ? tmp[i] & 15 : (tmp[i] & 31) + 9);
                    ++i;
                    --point;
                }
                if (i < (int)tmp.length())
                {
                    ++i;
                    while (i < (int)tmp.length())
                    {
                        number.push_back(tmp[i] <= '9' ? tmp[i] & 15 : (tmp[i] & 31) + 9);
                        ++i;
                    }
                }
                else
                {
                    point = 0;
                }
            }
            else if (pow(10, (int)log10(base)) == base)
            {
                int len = log10(base), v = 1, dot = tmp.find('.');
                if (dot == -1)
                {
                    point = 0;
                }
                else
                {
                    tmp.erase(dot, 1);
                    tmp = std::string((len - dot % len) % len, '0') + tmp;
                    point = -(dot / len + (dot % len != 0));
                }
                for (int i = 0; i < (int)tmp.length(); ++i)
                {
                    if (i % len == 0)
                    {
                        number.push_back((tmp[i] & 15) * (v = 1));
                    }
                    else
                    {
                        number.back() += (tmp[i] & 15) * (v *= 10);
                    }
                }
            }
            shrink_to_fit();
        }
        void clear() _GLIBCXX_NOEXCEPT
        {
            number.clear();
            mark = false;
            point = 0;
            up_per_time = 10;
        }
        void reset_point(const int &__Val) _GLIBCXX_NOEXCEPT
        {
            if (point < __Val)
            {
                *this >>= __Val - point;
            }
            if (point > __Val)
            {
                *this <<= point - __Val;
            }
            point = __Val;
        }
        void delete_left_zero() _GLIBCXX_NOEXCEPT
        {
            int __Val = 0;
            while (__Val < (int)size() && !number[__Val])
            {
                ++__Val;
            }
            reset_point(__Val + point);
        }
        void reset_up_per_time(const int &__Val) _GLIBCXX_NOEXCEPT
        {
            up_per_time = __Val;
        }
        bool empty() const _GLIBCXX_NOEXCEPT
        {
            return number.empty();
        }
        size_t size() const _GLIBCXX_NOEXCEPT
        {
            return number.size();
        }
        bool operator<(const bigfloat &__Val) const _GLIBCXX_NOEXCEPT
        {
            if (mark != __Val.mark)
            {
                return mark && !__Val.mark;
            }
            if (mark)
            {
                return -__Val < -*this;
            }
            if (empty() ^ __Val.empty())
            {
                return empty();
            }
            if (empty() && __Val.empty())
            {
                return false;
            }
            if (size() + point != __Val.size() + __Val.point)
            {
                return (int)size() + point < (int)__Val.size() + __Val.point;
            }
            for (int i = size() - 1, j = __Val.size() - 1; i >= 0 && j >= 0; --i, --j)
            {
                if (number[i] != __Val.number[j])
                {
                    return number[i] < __Val.number[j];
                }
            }
            if (size() != __Val.size())
            {
                return size() < __Val.size();
            }
            return false;
        }
        bool operator==(const bigfloat &__Val) const _GLIBCXX_NOEXCEPT
        {
            return mark == __Val.mark && point == __Val.point && number == __Val.number;
        }
        bool operator>(const bigfloat &__Val) const _GLIBCXX_NOEXCEPT
        {
            return __Val < *this;
        }
        bool operator<=(const bigfloat &__Val) const _GLIBCXX_NOEXCEPT
        {
            return *this == __Val || *this < __Val;
        }
        bool operator>=(const bigfloat &__Val) const _GLIBCXX_NOEXCEPT
        {
            return *this == __Val || *this > __Val;
        }
        bool operator!=(const bigfloat &__Val) const _GLIBCXX_NOEXCEPT
        {
            return !(*this == __Val);
        }
        bigfloat &operator++() _GLIBCXX_NOEXCEPT
        {
            if (empty())
            {
                number.push_back(1);
                return *this;
            }
            number.push_back(0);
            if (mark)
            {
                --number.front();
                for (int i = 0; number[i] < 0; ++i)
                {
                    number[i] = base - 1;
                    --number[i + 1];
                }
            }
            else
            {
                ++number.front();
                for (int i = 0; number[i] == base; ++i)
                {
                    number[i] = 0;
                    ++number[i + 1];
                }
            }
            shrink_to_fit();
            return *this;
        }
        bigfloat &operator--() _GLIBCXX_NOEXCEPT
        {
            if (empty())
            {
                number.push_back(1);
                mark = true;
                return *this;
            }
            number.push_back(0);
            if (mark)
            {
                ++number.front();
                for (int i = 0; number[i] == base; ++i)
                {
                    number[i] = 0;
                    ++number[i + 1];
                }
            }
            else
            {
                --number.front();
                for (int i = 0; number[i] < 0; ++i)
                {
                    number[i] = base - 1;
                    --number[i + 1];
                }
            }
            shrink_to_fit();
            return *this;
        }
        bigfloat operator++(int) _GLIBCXX_NOEXCEPT
        {
            bigfloat res = *this;
            ++*this;
            return res;
        }
        bigfloat operator--(int) _GLIBCXX_NOEXCEPT
        {
            bigfloat res = *this;
            --*this;
            return res;
        }
        bigfloat operator+() const _GLIBCXX_NOEXCEPT
        {
            return *this;
        }
        bigfloat operator-() const _GLIBCXX_NOEXCEPT
        {
            if (empty())
            {
                return *this;
            }
            bigfloat res = *this;
            res.mark ^= true;
            return res;
        }
        bigfloat operator+(const bigfloat &__Val) const _GLIBCXX_NOEXCEPT
        {
            if (empty())
            {
                return __Val;
            }
            if (__Val.empty())
            {
                return *this;
            }
            if (mark)
            {
                return __Val - (-*this);
            }
            if (__Val.mark)
            {
                return *this - (-__Val);
            }
            bigfloat res;
            res.number.assign(std::max(size() + point, __Val.size() - __Val.point) - std::min(point, __Val.point), 0);
            res.point = std::min(point, __Val.point);
            for (int i = 0; i < (int)size(); ++i)
            {
                res.number[i + point - res.point] += number[i];
            }
            for (int i = 0; i < (int)__Val.size(); ++i)
            {
                res.number[i + __Val.point - res.point] += __Val.number[i];
            }
            for (int i = 1; i < (int)res.size(); ++i)
            {
                if (res.number[i - 1] >= base)
                {
                    res.number[i - 1] -= base;
                    ++res.number[i];
                }
            }
            if (res.number.back() >= base)
            {
                res.number.back() -= base;
                res.number.push_back(1);
            }
            res.shrink_to_fit();
            return res;
        }
        bigfloat operator-(const bigfloat &__Val) const _GLIBCXX_NOEXCEPT
        {
            if (empty())
            {
                return -__Val;
            }
            if (__Val.empty())
            {
                return *this;
            }
            if (mark)
            {
                return -(-*this + __Val);
            }
            if (__Val.mark)
            {
                return *this + __Val;
            }
            if (*this < __Val)
            {
                return -(__Val - *this);
            }
            bigfloat res = *this;
            int move = 0;
            if (__Val.point < point)
            {
                res.reset_point(__Val.point);
            }
            else
            {
                move = __Val.point - res.point;
            }
            for (int i = 0; i < (int)__Val.size(); ++i)
            {
                if (res.number[i + move] < __Val.number[i])
                {
                    res.number[i + move] += base;
                    --res.number[i + move + 1];
                }
                res.number[i + move] -= __Val.number[i];
            }
            res.shrink_to_fit();
            return res;
        }
        bigfloat operator*(const bigfloat &__Val) const _GLIBCXX_NOEXCEPT
        {
            if (empty() || __Val.empty())
            {
                return 0;
            }
            bigfloat res;
            if (size() + __Val.size() <= 512)
            {
                res.mark = mark ^ __Val.mark;
                res.point = point + __Val.point;
                res.number.assign(size() + __Val.size(), 0);
                for (int i = 0; i < (int)size(); ++i)
                {
                    for (int j = 0; j < (int)__Val.size(); ++j)
                    {
                        res.number[i + j] += number[i] * __Val.number[j];
                    }
                }
                for (int i = 1; i < (int)res.size(); ++i)
                {
                    if (res.number[i - 1] >= base)
                    {
                        res.number[i] += res.number[i - 1] / base;
                        res.number[i - 1] %= base;
                    }
                }
                res.shrink_to_fit();
            }
            else if (size() + __Val.size() <= 40000000)
            {
                res.point = point + __Val.point;
                res.mark = mark ^ __Val.mark;
                unsigned len = 32 - __builtin_clz(size() + __Val.size()), siz = 1 << len;
                std::vector<complex<double> /**/> A(siz), B;
                for (unsigned i = 0; i != size(); ++i)
                    A[i].x = number[i];
                for (unsigned i = 0; i != __Val.size(); ++i)
                    A[i].y = __Val.number[i];
                static std::vector<int> rev;
                if (rev.size() != siz)
                {
                    rev.resize(siz);
                    for (unsigned i = 0; i != siz; ++i)
                        rev[i] = rev[i >> 1] | ((i & 1) << len);
                }
                for (unsigned i = 0; i != siz; ++i)
                    if (i < rev[i])
                        std::swap(A[i], A[rev[i]]);
                for (unsigned Len = 1; Len != siz; Len = Len + Len)
                {
                    static complex<double> Y1, Y2;
                    complex<double> w1(cos(PI / Len), sin(PI / Len));
                    for (unsigned i = 0; i != len; i += Len + Len)
                    {
                        complex<double> w(1, 0);
                        for (unsigned j = 0; j != Len; ++j, w *= w1)
                        {
                            Y1 = A[i + j], Y2 = w * A[i + j + Len];
                            A[i + j /* */] = Y1 + Y2;
                            A[i + j + Len] = Y1 - Y2;
                        }
                    }
                }
                for (std::vector<complex<double> /**/>::iterator i = A.begin(); i != A.end(); ++i)
                    *i *= *i;
                for (unsigned i = 0; i != A.size(); ++i)
                    B.push_back(A[-i & (siz - 1)] - conj(A[i]));
                for (unsigned i = 0; i != siz; ++i)
                    if (i < rev[i])
                        std::swap(B[i], B[rev[i]]);
                for (unsigned Len = 1; Len != siz; Len = Len + Len)
                {
                    static complex<double> Y1, Y2;
                    complex<double> w1(cos(PI / Len), sin(PI / Len));
                    for (unsigned i = 0; i != len; i += Len + Len)
                    {
                        complex<double> w(1, 0);
                        for (unsigned j = 0; j != Len; ++j, w *= w1)
                        {
                            Y1 = B[i + j], Y2 = w * B[i + j + Len];
                            B[i + j /* */] = Y1 + Y2;
                            B[i + j + Len] = Y1 - Y2;
                        }
                    }
                }
                for (int i = 0; i < (int)B.size(); ++i)
                {
                    res.number.push_back(B[i].y / (siz << 2) + 0.5);
                }
                for (int i = 1; i < (int)res.size(); ++i)
                {
                    if (res.number[i - 1] >= base)
                    {
                        res.number[i] += res.number[i - 1] / base;
                        res.number[i - 1] %= base;
                    }
                }
                if (res.number.back() >= base)
                {
                    int up = res.number.back() / base;
                    res.number.back() %= base;
                    res.number.push_back(up);
                }
                res.shrink_to_fit();
            }
            return res;
        }
        bigfloat inv(const int &accuracy = -100) const _GLIBCXX_NOEXCEPT
        {
            const double start_waste = log2(size()) + 0.55;
            __glibcxx_requires_subscript(*this != 0);
            bigfloat res;
            res.point = -(point + (int)size());
            res.mark = mark;
            res.number.push_back(0);
            for (int i = 30; i >= 0; --i)
            {
                if (number.back() * (res.number.front() | 1 << i) < (long long)base)
                {
                    res.number.front() |= 1 << i;
                }
            }
            double waste = start_waste;
            while (res.point + 2000 > accuracy)
            {
                res = res * (-*this * res + 2);
                waste *= 2;
                if (waste > 100)
                {
                    res.reset_point(res.point + (int)waste);
                    waste = start_waste;
                }
            }
            res.reset_point(accuracy);
            res = res * (-*this * res + 2);
            res.reset_point(accuracy);
            return res;
        }
        bigfloat operator/(const bigfloat &__Val) const _GLIBCXX_NOEXCEPT
        {
            return *this * __Val.inv(-((int)size() - (int)__Val.size() + up_per_time));
        }
        bigfloat operator/(const long long &__Val) const _GLIBCXX_NOEXCEPT
        {
            __glibcxx_requires_subscript(__Val);
            if (empty())
            {
                return *this;
            }
            bigfloat res;
            long long v = __Val < 0 ? -__Val : __Val, Left = 0;
            res.mark = mark ^ (__Val < 0);
            res.point = point - up_per_time;
            res.number.assign(size() + up_per_time, 0);
            for (int i = size() - 1; i >= 0; --i)
            {
                Left = Left * base + number[i];
                res.number[i + up_per_time] = Left / v;
                Left %= v;
            }
            for (int i = -1; i >= -up_per_time; --i)
            {
                Left = Left * base;
                res.number[i + up_per_time] = Left / v;
                Left %= v;
            }
            res.shrink_to_fit();
            return res;
        }
        bigfloat operator<<(const int &__Val) const _GLIBCXX_NOEXCEPT
        {
            bigfloat res;
            res.mark = mark;
            res.point = point;
            res.number.resize(size() + __Val);
            std::copy(number.begin(), number.end(), res.number.begin() + __Val);
            return res;
        }
        bigfloat operator>>(const int &__Val) const _GLIBCXX_NOEXCEPT
        {
            if (__Val >= (int)size())
            {
                return 0;
            }
            bigfloat res;
            res.mark = mark;
            res.point = point;
            res.number.resize(size() - __Val);
            std::copy(number.begin() + __Val, number.end(), res.number.begin());
            return res;
        }
        friend bigfloat abs(const bigfloat &__Val) _GLIBCXX_NOEXCEPT
        {
            bigfloat res = __Val;
            res.mark = false;
            return res;
        }
        friend std::istream &operator>>(std::istream &ist, bigfloat &__Val) _GLIBCXX_NOEXCEPT
        {
            __Val.clear();
            std::string tmp;
            ist >> tmp;
            __Val = tmp;
            return ist;
        }
        friend std::ostream &operator<<(std::ostream &ost, const bigfloat &__Val)
        {
            return ost << std::string(__Val);
        }
        operator std::string() const
        {
            std::string res;
            if (empty())
            {
                return "0";
            }
            if (mark)
            {
                res = "-";
            }
            bool dot = false;
            if (base <= 36)
            {
                if ((int)size() + point <= 0)
                {
                    res += "0.";
                    dot = true;
                    for (int i = 0; i > (int)size() + point; --i)
                        res += '0';
                }
                for (int i = size() - 1; i >= 0; --i)
                {
                    res += number[i] < 10 ? number[i] | 48 : 64 | (number[i] - 9);
                    if (!(i + point) && i)
                    {
                        res += '.';
                        dot = true;
                    }
                }
                for (int i = point; i > 0; --i)
                    res += '0';
            }
            else if (power(10, (int)log10(base)) == base)
            {
                if ((int)size() + point <= 0)
                {
                    res += "0.";
                    dot = true;
                    for (int i = 0; i > (int)size() + point; --i)
                        for (int j = base / 10; j > 0; j /= 10)
                            res += '0';
                }
                for (int i = base / 10; i; i /= 10)
                    if (dot || number.back() >= i)
                        res += '0' + number.back() / i % 10;
                for (int i = (int)size() - 2; i >= 0; --i)
                {
                    if (i + point == -1)
                    {
                        res += '.';
                        dot = true;
                    }
                    for (int j = base / 10; j; j /= 10)
                        res += '0' + number[i] / j % 10;
                }
            }
            while (dot && *--res.end() == '0')
                res.erase(--res.end());
            if (*--res.end() == '.')
                res.erase(--res.end());
            return res;
        }
        bigfloat &operator+=(const bigfloat &__Val) _GLIBCXX_NOEXCEPT { return *this = *this + __Val; }
        bigfloat &operator-=(const bigfloat &__Val) _GLIBCXX_NOEXCEPT { return *this = *this - __Val; }
        bigfloat &operator*=(const bigfloat &__Val) _GLIBCXX_NOEXCEPT { return *this = *this * __Val; }
        bigfloat &operator/=(const bigfloat &__Val) _GLIBCXX_NOEXCEPT { return *this = *this / __Val; }
        bigfloat &operator/=(const long long &__Val) _GLIBCXX_NOEXCEPT { return *this = *this / __Val; }
        bigfloat &operator<<=(const int &__Val) _GLIBCXX_NOEXCEPT { return *this = *this << __Val; }
        bigfloat &operator>>=(const int &__Val) _GLIBCXX_NOEXCEPT { return *this = *this >> __Val; }
    };
    template <size_t base = 10>
    class bigint
    {
    private:
        std::vector<long long> number;
        bool mark;
        void shrink_to_fit() _GLIBCXX_NOEXCEPT
        {
            while (!(number.empty() || number.back()))
                number.pop_back();
            if (number.empty())
                mark = false;
        }

    public:
        class iterator
        {
        private:
            std::vector<long long>::iterator Iter;

        public:
            _GLIBCXX14_CONSTEXPR iterator() _GLIBCXX_NOEXCEPT {}
            iterator(const std::vector<long long>::iterator &__Val) _GLIBCXX_NOEXCEPT { Iter = __Val; }
            iterator operator++() _GLIBCXX_NOEXCEPT { return ++Iter; }
            iterator operator++(int) _GLIBCXX_NOEXCEPT { return Iter++; }
            iterator operator--() _GLIBCXX_NOEXCEPT { return --Iter; }
            iterator operator--(int) _GLIBCXX_NOEXCEPT { return Iter--; }
            bool operator==(const iterator &__Val) const _GLIBCXX_NOEXCEPT { return Iter == __Val.Iter; }
            bool operator!=(const iterator &__Val) const _GLIBCXX_NOEXCEPT { return Iter != __Val.Iter; }
            long long &operator*() const _GLIBCXX_NOEXCEPT { return *Iter; }
        };
        class reverse_iterator
        {
        private:
            std::vector<long long>::reverse_iterator Iter;

        public:
            _GLIBCXX14_CONSTEXPR reverse_iterator() _GLIBCXX_NOEXCEPT {}
            reverse_iterator(const std::vector<long long>::reverse_iterator &__Val) _GLIBCXX_NOEXCEPT { Iter = __Val; }
            reverse_iterator operator++() _GLIBCXX_NOEXCEPT { return ++Iter; }
            reverse_iterator operator++(int) _GLIBCXX_NOEXCEPT { return Iter++; }
            reverse_iterator operator--() _GLIBCXX_NOEXCEPT { return --Iter; }
            reverse_iterator operator--(int) _GLIBCXX_NOEXCEPT { return Iter--; }
            bool operator==(const reverse_iterator &__Val) const _GLIBCXX_NOEXCEPT { return Iter == __Val.Iter; }
            bool operator!=(const reverse_iterator &__Val) const _GLIBCXX_NOEXCEPT { return Iter != __Val.Iter; }
            long long &operator*() const _GLIBCXX_NOEXCEPT { return *Iter; }
        };
        iterator begin() _GLIBCXX_NOEXCEPT
        {
            return number.begin();
        }
        iterator end() _GLIBCXX_NOEXCEPT
        {
            return number.end();
        }
        reverse_iterator rbegin() _GLIBCXX_NOEXCEPT
        {
            return number.rbegin();
        }
        reverse_iterator rend() _GLIBCXX_NOEXCEPT
        {
            return number.rend();
        }
        _GLIBCXX14_CONSTEXPR bigint() _GLIBCXX_NOEXCEPT : number(), mark() {}
        template <typename __Tp>
        bigint(const __Tp &__Val) _GLIBCXX_NOEXCEPT : number(), mark()
        {
            __Tp tmp = __Val;
            if (tmp < 0)
            {
                mark = true;
                tmp = -tmp;
            }
            while (tmp > 0)
            {
                number.push_back(tmp % base);
                tmp /= base;
            }
        }
        bigint(const std::string &__Val) _GLIBCXX_NOEXCEPT : number(), mark()
        {
            std::string tmp = __Val;
            std::reverse(tmp.begin(), tmp.end());
            if (*--tmp.end() == '-')
            {
                mark = true;
                tmp.erase(--tmp.end());
            }
            if (base <= 36)
                for (int i = 0; i < (int)tmp.length(); ++i)
                    number.push_back(tmp[i] <= '9' ? tmp[i] & 15 : (tmp[i] & 31) + 9);
            else if (pow(10, (int)log10(base)) == base)
                for (int i = 0, len = log10(base), v = 1; i < (int)tmp.length(); ++i)
                {
                    if (i % len == 0)
                        number.push_back((tmp[i] & 15) * (v = 1));
                    else
                        number.back() += (tmp[i] & 15) * (v *= 10);
                }
            shrink_to_fit();
        }
        void clear() _GLIBCXX_NOEXCEPT
        {
            number.clear();
            mark = false;
        }
        bool empty() const _GLIBCXX_NOEXCEPT { return number.empty(); }
        bool operator<(const bigint &__Val) const _GLIBCXX_NOEXCEPT
        {
            if (mark != __Val.mark)
                return mark && !__Val.mark;
            if (mark)
                return -__Val < -*this;
            if (empty() ^ __Val.empty())
                return empty();
            if (empty() && __Val.empty())
                return false;
            if (size() != __Val.size())
                return size() < __Val.size();
            for (int i = size() - 1; i >= 0; --i)
                if (number[i] != __Val.number[i])
                    return number[i] < __Val.number[i];
            return false;
        }
        bool operator!=(const bigint &__Val) const _GLIBCXX_NOEXCEPT { return (*this < __Val) || (__Val < *this); }
        bool operator==(const bigint &__Val) const _GLIBCXX_NOEXCEPT { return !((*this < __Val) || (__Val < *this)); }
        bool operator>(const bigint &__Val) const _GLIBCXX_NOEXCEPT { return __Val < *this; }
        bool operator<=(const bigint &__Val) const _GLIBCXX_NOEXCEPT { return !(__Val < *this); }
        bool operator>=(const bigint &__Val) const _GLIBCXX_NOEXCEPT { return !(*this < __Val); }
        size_t size() const _GLIBCXX_NOEXCEPT { return number.size(); }
        bigint &operator++() _GLIBCXX_NOEXCEPT
        {
            if (empty())
            {
                number.push_back(1);
                return *this;
            }
            number.push_back(0);
            if (mark)
            {
                --number.front();
                for (int i = 0; number[i] < 0; ++i)
                {
                    number[i] = base - 1;
                    --number[i + 1];
                }
            }
            else
            {
                ++number.front();
                for (int i = 0; number[i] == base; ++i)
                {
                    number[i] = 0;
                    ++number[i + 1];
                }
            }
            shrink_to_fit();
            return *this;
        }
        bigint &operator--() _GLIBCXX_NOEXCEPT
        {
            if (empty())
            {
                number.push_back(1);
                mark = true;
                return *this;
            }
            number.push_back(0);
            if (mark)
            {
                ++number.front();
                for (int i = 0; number[i] == base; ++i)
                {
                    number[i] = 0;
                    ++number[i + 1];
                }
            }
            else
            {
                --number.front();
                for (int i = 0; number[i] < 0; ++i)
                {
                    number[i] = base - 1;
                    --number[i + 1];
                }
            }
            shrink_to_fit();
        }
        bigint operator++(int) _GLIBCXX_NOEXCEPT
        {
            bigint res = *this;
            ++*this;
            return res;
        }
        bigint operator--(int) _GLIBCXX_NOEXCEPT
        {
            bigint res = *this;
            --*this;
            return res;
        }
        bigint operator+() const _GLIBCXX_NOEXCEPT { return *this; }
        bigint operator-() const _GLIBCXX_NOEXCEPT
        {
            if (empty())
            {
                return *this;
            }
            bigint res = *this;
            res.mark ^= true;
            return res;
        }
        bigint operator+(const bigint &__Val) const _GLIBCXX_NOEXCEPT
        {
            if (empty())
                return __Val;
            if (__Val.empty())
                return *this;
            if (mark)
                return __Val - (-*this);
            if (__Val.mark)
                return *this - (-__Val);
            bigint res;
            res.number.assign(std::max(size(), __Val.size()), 0);
            for (int i = 0; i < (int)size(); ++i)
                res.number[i] += number[i];
            for (int i = 0; i < (int)__Val.size(); ++i)
                res.number[i] += __Val.number[i];
            for (int i = 1; i < (int)res.size(); ++i)
                if (res.number[i - 1] >= base)
                {
                    res.number[i - 1] -= base;
                    ++res.number[i];
                }
            if (res.number.back() >= base)
            {
                res.number.back() -= base;
                res.number.push_back(1);
            }
            return res;
        }
        bigint operator-(const bigint &__Val) const _GLIBCXX_NOEXCEPT
        {
            if (empty())
                return -__Val;
            if (__Val.empty())
                return *this;
            if (mark)
                return -(-*this + __Val);
            if (__Val.mark)
                return *this + __Val;
            if (*this < __Val)
                return -(__Val - *this);
            bigint res = *this;
            for (int i = 0; i < (int)__Val.size(); ++i)
            {
                if (res.number[i] < __Val.number[i])
                {
                    res.number[i] += base;
                    --res.number[i + 1];
                }
                res.number[i] -= __Val.number[i];
            }
            for (int i = __Val.size(); res.number[i] < 0; ++i)
            {
                res.number[i] += base;
                --res.number[i + 1];
            }
            res.shrink_to_fit();
            return res;
        }
        bigint operator*(const bigint &__Val) const _GLIBCXX_NOEXCEPT
        {
            if (empty() || __Val.empty())
                return 0;
            bigint res;
            if (size() + __Val.size() <= 512)
            {
                res.mark = mark ^ __Val.mark;
                res.number.assign(size() + __Val.size(), 0);
                for (int i = 0; i < (int)size(); ++i)
                    for (int j = 0; j < (int)__Val.size(); ++j)
                        res.number[i + j] += number[i] * __Val.number[j];
                for (int i = 1; i < (int)res.size(); ++i)
                    if (res.number[i - 1] >= base)
                    {
                        res.number[i] += res.number[i - 1] / base;
                        res.number[i - 1] %= base;
                    }
                res.shrink_to_fit();
            }
            else if (size() + __Val.size() <= 100000)
            {
                int len = std::min(size(), __Val.size()) >> 1;
                bigint AC = (*this >> len) * (__Val >> len);
                bigint BD = (*this & len) * (__Val & len);
                res = (AC << (len << 1)) +
                      ((((*this >> len) + (*this & len)) * ((__Val >> len) + (__Val & len)) - AC - BD) << len) +
                      BD;
            }
            else if (size() + __Val.size() <= 40000000)
            {
                res.mark = mark ^ __Val.mark;
                unsigned len = 32 - __builtin_clz(size() + __Val.size()), siz = 1 << len;
                std::vector<complex<double> /**/> A(siz), B;
                for (unsigned i = 0; i != size(); ++i)
                    A[i].x = number[i];
                for (unsigned i = 0; i != __Val.size(); ++i)
                    A[i].y = __Val.number[i];
                static std::vector<int> rev;
                if (rev.size() != siz)
                {
                    rev.resize(siz);
                    for (unsigned i = 0; i != siz; ++i)
                        rev[i] = rev[i >> 1] | ((i & 1) << len);
                }
                for (unsigned i = 0; i != siz; ++i)
                    if (i < rev[i])
                        std::swap(A[i], A[rev[i]]);
                for (unsigned Len = 1; Len != siz; Len = Len + Len)
                {
                    static complex<double> Y1, Y2;
                    complex<double> w1(cos(PI / Len), sin(PI / Len));
                    for (unsigned i = 0; i != len; i += Len + Len)
                    {
                        complex<double> w(1.0, 0.0);
                        for (unsigned j = 0; j != Len; ++j, w *= w1)
                        {
                            Y1 = A[i + j], Y2 = w * A[i + j + Len];
                            A[i + j /* */] = Y1 + Y2;
                            A[i + j + Len] = Y1 - Y2;
                        }
                    }
                }
                for (std::vector<complex<double> /**/>::iterator i = A.begin(); i != A.end(); ++i)
                    *i *= *i;
                for (unsigned i = 0; i != A.size(); ++i)
                    B.push_back(A[-i & (siz - 1)] - conj(A[i]));
                for (unsigned i = 0; i != siz; ++i)
                    if (i < rev[i])
                        std::swap(B[i], B[rev[i]]);
                for (unsigned Len = 1; Len != siz; Len = Len + Len)
                {
                    static complex<double> Y1, Y2;
                    complex<double> w1(cos(PI / Len), sin(PI / Len));
                    for (unsigned i = 0; i != len; i += Len + Len)
                    {
                        complex<double> w(1.0, 0.0);
                        for (unsigned j = 0; j != Len; ++j, w *= w1)
                        {
                            Y1 = B[i + j], Y2 = w * B[i + j + Len];
                            B[i + j /* */] = Y1 + Y2;
                            B[i + j + Len] = Y1 - Y2;
                        }
                    }
                }
                for (int i = 0; i < (int)B.size(); ++i)
                {
                    res.number.push_back(B[i].y / (siz << 2) + 0.5);
                }
                for (int i = 1; i < (int)res.size(); ++i)
                {
                    if (res.number[i - 1] >= base)
                    {
                        res.number[i] += res.number[i - 1] / base;
                        res.number[i - 1] %= base;
                    }
                }
                if (res.number.back() >= base)
                {
                    int up = res.number.back() / base;
                    res.number.back() %= base;
                    res.number.push_back(up);
                }
                res.shrink_to_fit();
            }
            return res;
        }
        bigint operator/(const bigint &__Val) const _GLIBCXX_NOEXCEPT
        {
            __glibcxx_requires_subscript(__Val != 0);
            if (mark)
            {
                return -(-*this / __Val);
            }
            if (__Val.mark)
            {
                return -(*this / -__Val);
            }
            bigint res = bigfloat<base>(*this) / bigfloat<base>(__Val);
            if (__Val * (res + 1) <= *this)
            {
                res += 1;
            }
            return res;
        }
        bigint operator/(const long long &__Val) const _GLIBCXX_NOEXCEPT
        {
            __glibcxx_requires_subscript(__Val);
            if (empty())
            {
                return *this;
            }
            bigint res;
            long long v = __Val < 0 ? -__Val : __Val, Left = 0;
            res.mark = mark ^ (__Val < 0);
            res.number.assign(size(), 0);
            for (int i = res.size() - 1; i >= 0; --i)
            {
                Left = Left * base + number[i];
                res.number[i] = Left / v;
                Left %= v;
            }
            res.shrink_to_fit();
            return res;
        }
        long long operator%(const long long &__Val) const _GLIBCXX_NOEXCEPT
        {
            __glibcxx_requires_subscript(__Val);
            long long Left = 0;
            for (int i = size() - 1; i >= 0; --i)
            {
                Left = Left * base + number[i];
                Left %= __Val;
            }
            return Left * (mark ? -1 : 1);
        }
        bigint operator<<(const int &__Val) const _GLIBCXX_NOEXCEPT
        {
            bigint res;
            res.mark = mark;
            res.number.resize(size() + __Val);
            std::copy(number.begin(), number.end(), res.number.begin() + __Val);
            return res;
        }
        bigint operator>>(const int &__Val) const _GLIBCXX_NOEXCEPT
        {
            if (__Val >= (int)size())
            {
                return 0;
            }
            bigint res;
            res.mark = mark;
            res.number.resize(size() - __Val);
            std::copy(number.begin() + __Val, number.end(), res.number.begin());
            return res;
        }
        bigint operator&(const int &__Val) const _GLIBCXX_NOEXCEPT
        {
            bigint res = *this;
            res.number.resize(__Val);
            return res;
        }
        bigint operator|(const bigint &__Val) const _GLIBCXX_NOEXCEPT
        {
            bigint res = __Val;
            res.mark |= mark;
            res.number.resize(size() + __Val.size());
            std::copy(number.begin(), number.end(), res.number.begin() + __Val.size());
            return res;
        }
        friend bigint abs(const bigint &__Val) _GLIBCXX_NOEXCEPT
        {
            bigint res = __Val;
            res.mark = false;
            return res;
        }
        friend std::istream &operator>>(std::istream &ist, bigint &__Val) _GLIBCXX_NOEXCEPT
        {
            __Val.clear();
            std::string tmp;
            ist >> tmp;
            __Val = tmp;
            return ist;
        }
        friend std::ostream &operator<<(std::ostream &ost, const bigint &__Val) _GLIBCXX_NOEXCEPT
        {
            return ost << std::string(__Val);
        }
        operator std::string() const
        {
            if (empty())
                return "0";
            std::string res;
            if (mark)
                res = "-";
            if (base <= 36)
                for (int i = size() - 1; i >= 0; --i)
                    res += (char)(number[i] < 10 ? number[i] | 48 : 64 | (number[i] - 9));
            else if (pow(10, (int)log10(base)) == base)
            {
                for (int i = power(10, (int)log10(number.back())); i; i /= 10)
                    res += '0' + number.back() / i % 10;
                for (int i = (int)size() - 2; i >= 0; --i)
                    for (int j = base / 10; j > number[i]; j /= 10)
                        res += '0' + number[i] / j % 10;
            }
            return res;
        }
        bigint &operator+=(const bigint &__Val) _GLIBCXX_NOEXCEPT
        {
            return *this = *this + __Val;
        }
        bigint &operator-=(const bigint &__Val) _GLIBCXX_NOEXCEPT
        {
            return *this = *this - __Val;
        }
        bigint &operator*=(const bigint &__Val) _GLIBCXX_NOEXCEPT
        {
            return *this = *this * __Val;
        }
        bigint &operator*=(const long long &__Val) _GLIBCXX_NOEXCEPT
        {
            return *this = *this * __Val;
        }
        bigint &operator/=(const bigint &__Val) _GLIBCXX_NOEXCEPT
        {
            return *this = *this / __Val;
        }
        bigint &operator/=(const long long &__Val) _GLIBCXX_NOEXCEPT
        {
            return *this = *this / __Val;
        }
        bigint &operator<<=(const int &__Val) _GLIBCXX_NOEXCEPT
        {
            return *this = *this << __Val;
        }
        bigint &operator>>=(const int &__Val) _GLIBCXX_NOEXCEPT
        {
            return *this = *this >> __Val;
        }
    };
    template <typename __Tp>
    __Tp power(__Tp A, long long B)
    {
        __Tp res = 1;
        while (B)
        {
            if (B & 1)
            {
                res *= A;
            }
            A *= A;
            B >>= 1;
        }
        return res;
    }
} // namespace zjcSTL