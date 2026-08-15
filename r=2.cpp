/**
 * 序列搜索器 r=2 —— 基础版（无树状日志）
 * 
 * 功能：
 * - 深度优先回溯搜索严格递增正整数序列，满足 T_{n+1} = 2*T_n
 * - 应用多种数学约束剪枝
 * - 迭代回溯填充新增段，避免递归栈溢出
 * - 日志系统解耦，支持详细文件日志、统计日志、空日志
 * - 保存完整 promising 序列快照，控制台输出完整前缀
 * 
 * 日志模式切换：
 *   - 详细文件日志：init_detailed_log("r2.log");
 *   - 统计日志：     current_sink = &stats_sink;
 *   - 空日志：       current_sink = &null_sink;
 * 
 * 编译：需要 GCC/Clang（支持 __builtin_add_overflow / __builtin_mul_overflow）
 *      链接数学库 -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <stddef.h>

#ifdef _WIN32
    #include <windows.h>
    #define SLEEP_MS(ms) Sleep(ms)
#else
    #include <unistd.h>
    #define SLEEP_MS(ms) usleep((ms) * 1000)
#endif

typedef long long ll;

#define LIMIT (1LL << 60)                  // T_n 超过此值视为严格 promising
#define MAX_VAL 1000000000LL               // 数值安全上限，防止 long long 乘法溢出
#define PROMISING_N_THRESHOLD 20           // 扩展步数启发式阈值
#define PROMISING_ITEM_THRESHOLD 100000000LL // 末项启发式阈值

/* ---------------- 失败原因枚举 ---------------- */
/*
 * 失败原因枚举：每条判定分支都会记录一次失败，便于日志统计与剪枝调试。
 * 设计目标是把“为什么这条分支不可能扩展”拆成可观察的几类条件，
 * 这样后续修改算法时可以直接定位到对应的判定逻辑。
 */
typedef enum {
    OK = 0,
    FAIL_STEP_BELOW_MIN,
    FAIL_STEP_ABOVE_MAX,
    FAIL_REFINED_INEQ,
    FAIL_MIN_SUM_EXCEED,
    FAIL_FILL_IMPOSSIBLE,
    FAIL_L_EXTRA,
    FAIL_NEED_SUM_NOT_ZERO,
    FAIL_INITIAL_RANGE_EMPTY,
    FAIL_LIMIT_REACHED
} FailureReason;

const char* failure_reason_str(FailureReason r) {
    switch (r) {
        case OK: return "OK";
        case FAIL_STEP_BELOW_MIN: return "步长小于2";
        case FAIL_STEP_ABOVE_MAX: return "步长超过上界";
        case FAIL_REFINED_INEQ: return "精细不等式不成立";
        case FAIL_MIN_SUM_EXCEED: return "新增段最小和超过Tn";
        case FAIL_FILL_IMPOSSIBLE: return "填充范围为空";
        case FAIL_L_EXTRA: return "L取上界额外限制失败";
        case FAIL_NEED_SUM_NOT_ZERO: return "填充后need_sum非零";
        case FAIL_INITIAL_RANGE_EMPTY: return "初始段枚举范围为空";
        case FAIL_LIMIT_REACHED: return "达到LIMIT(promising)";
        default: return "未知";
    }
}

/* ---------------- 状态结构 ---------------- */
/*
 * State 保存当前搜索前缀的全部必要信息。
 * a[1..m] 表示已经确定好的序列前缀，m 是当前最后一项下标。
 * n 表示当前已确认的长度（严格来说是满足 T_n 的项数），
 * Tn 表示当前前缀对应的和，也就是 T_n。
 *
 * 整体搜索思路：
 *   1. 先枚举初始段（a1 与 a2 之后的部分）；
 *   2. 通过 try_extend 试图在末尾继续扩展新的区间；
 *   3. 若满足某些启发式条件，则保存 promising 序列快照并提前结束。
 */
typedef struct {
    ll *a;            // 动态数组，索引从 1 开始
    ll m;             // 当前已填充的最大索引（即 a_n）
    ll n;             // 当前步数（已满足 T_n 的个数）
    ll Tn;            // 当前和 T_n
    size_t cap;       // 数组容量
    ll N_threshold;   // 简化全局上界阈值（未使用，保留）
} State;

State st;

/* 保存 promising 序列完整快照 */
ll *promising_a = NULL;
ll promising_m = 0;

/* ---------------- 日志系统接口 ---------------- */
typedef struct LogSink {
    void (*write)(FailureReason reason, const State *s, const char *where);
    void (*flush)(void);
    void (*close)(void);
} LogSink;

extern const LogSink null_sink;
const LogSink *current_sink = &null_sink;

/* ---------------- 工具函数 ---------------- */
static inline bool add_overflow(ll a, ll b, ll *res) {
    return __builtin_add_overflow(a, b, res);
}
static inline bool mul_overflow(ll a, ll b, ll *res) {
    return __builtin_mul_overflow(a, b, res);
}

void *safe_realloc(void *ptr, size_t size, State *s) {
    int retries = 0;
    const int max_retries = 10;
    void *new_ptr = NULL;
    while (retries < max_retries) {
        new_ptr = realloc(ptr, size);
        if (new_ptr) return new_ptr;
        retries++;
        fprintf(stderr, "内存分配失败 (第 %d 次)，暂停 5 秒后重试...\n", retries);
        SLEEP_MS(5000);
    }
    fprintf(stderr, "内存分配失败，无法继续。当前状态:\n");
    if (s->a != NULL) {
        fprintf(stderr, "  a1=%lld, n=%lld, m=%lld, Tn=%lld\n",
                s->a[1], s->n, s->m, s->Tn);
        if (s->m > 0) {
            fprintf(stderr, "  序列前缀: ");
            for (ll i = 1; i <= s->m && i <= 20; i++) {
                fprintf(stderr, "%lld ", s->a[i]);
            }
            fprintf(stderr, "\n");
        }
    } else {
        fprintf(stderr, "  (数组尚未分配)\n");
    }
    current_sink->flush();
    current_sink->close();
    exit(1);
}

void ensure_capacity(State *s, ll need) {
    if (need <= (ll)s->cap) return;
    size_t new_cap = s->cap ? s->cap * 2 : 16;
    while (new_cap < (size_t)need) new_cap *= 2;
    s->a = (ll *)safe_realloc(s->a, (new_cap + 1) * sizeof(ll), s);
    s->cap = new_cap;
}

void print_sequence_values(const State *s, FILE *out, ll max_items) {
    ll limit = (s->m < max_items) ? s->m : max_items;
    for (ll i = 1; i <= limit; i++) {
        fprintf(out, "%lld", s->a[i]);
        if (i < limit) fprintf(out, " ");
    }
    if (s->m > max_items) {
        fprintf(out, " ... (共 %lld 项)", s->m);
    }
    fprintf(out, "\n");
}

void save_promising_sequence(const State *s) {
    if (promising_a) {
        free(promising_a);
        promising_a = NULL;
    }
    promising_m = s->m;
    promising_a = (ll *)malloc((promising_m + 1) * sizeof(ll));
    if (!promising_a) {
        fprintf(stderr, "错误：无法分配 promising 序列内存\n");
        exit(1);
    }
    for (ll i = 1; i <= promising_m; i++) {
        promising_a[i] = s->a[i];
    }
}

/* ---------------- 日志函数 ---------------- */
void log_attempt(State *s, FailureReason reason, const char *where) {
    current_sink->write(reason, s, where);
}

/* ---------------- 空日志实现 ---------------- */
static void null_write(FailureReason reason, const State *s, const char *where) {
    (void)reason; (void)s; (void)where;
}
static void null_flush(void) {}
static void null_close(void) {}

const LogSink null_sink = {
    .write = null_write,
    .flush = null_flush,
    .close = null_close
};

/* ---------------- 详细文件日志实现 ---------------- */
static FILE *detailed_file = NULL;

static void detailed_write(FailureReason reason, const State *s, const char *where) {
    if (!detailed_file) return;
    fprintf(detailed_file, "[%s] a1=%lld, n=%lld, m=%lld, Tn=%lld  位置: %s\n",
            failure_reason_str(reason), s->a[1], s->n, s->m, s->Tn, where);
    fprintf(detailed_file, "  前缀: ");
    print_sequence_values(s, detailed_file, 20);
}
static void detailed_flush(void) {
    if (detailed_file) fflush(detailed_file);
}
static void detailed_close(void) {
    if (detailed_file) {
        fclose(detailed_file);
        detailed_file = NULL;
    }
}

LogSink detailed_sink = {
    .write = detailed_write,
    .flush = detailed_flush,
    .close = detailed_close
};

bool init_detailed_log(const char *path) {
    detailed_file = fopen(path, "w");
    if (!detailed_file) return false;
    setvbuf(detailed_file, NULL, _IOFBF, 4 << 20);
    current_sink = &detailed_sink;
    return true;
}

/* ---------------- 统计日志实现 ---------------- */
typedef struct {
    long long fail_step_below_min;
    long long fail_step_above_max;
    long long fail_refined_ineq;
    long long fail_min_sum_exceed;
    long long fail_fill_impossible;
    long long fail_L_extra;
    long long fail_need_sum_not_zero;
    long long fail_initial_range_empty;
    long long limit_reached;
} FailureStats;

FailureStats stats = {0};

static void stats_write(FailureReason reason, const State *s, const char *where) {
    (void)s; (void)where;
    switch (reason) {
        case FAIL_STEP_BELOW_MIN: stats.fail_step_below_min++; break;
        case FAIL_STEP_ABOVE_MAX: stats.fail_step_above_max++; break;
        case FAIL_REFINED_INEQ: stats.fail_refined_ineq++; break;
        case FAIL_MIN_SUM_EXCEED: stats.fail_min_sum_exceed++; break;
        case FAIL_FILL_IMPOSSIBLE: stats.fail_fill_impossible++; break;
        case FAIL_L_EXTRA: stats.fail_L_extra++; break;
        case FAIL_NEED_SUM_NOT_ZERO: stats.fail_need_sum_not_zero++; break;
        case FAIL_INITIAL_RANGE_EMPTY: stats.fail_initial_range_empty++; break;
        case FAIL_LIMIT_REACHED: stats.limit_reached++; break;
        default: break;
    }
}
static void stats_flush(void) {}
static void stats_close(void) {
    fprintf(stderr, "\n===== 失败统计 =====\n");
    fprintf(stderr, "步长小于2: %lld\n", stats.fail_step_below_min);
    fprintf(stderr, "步长超过上界: %lld\n", stats.fail_step_above_max);
    fprintf(stderr, "精细不等式不成立: %lld\n", stats.fail_refined_ineq);
    fprintf(stderr, "新增段最小和超过Tn: %lld\n", stats.fail_min_sum_exceed);
    fprintf(stderr, "填充范围为空: %lld\n", stats.fail_fill_impossible);
    fprintf(stderr, "L取上界额外限制失败: %lld\n", stats.fail_L_extra);
    fprintf(stderr, "填充后need_sum非零: %lld\n", stats.fail_need_sum_not_zero);
    fprintf(stderr, "初始段枚举范围为空: %lld\n", stats.fail_initial_range_empty);
    fprintf(stderr, "达到LIMIT(promising): %lld\n", stats.limit_reached);
}

LogSink stats_sink = {
    .write = stats_write,
    .flush = stats_flush,
    .close = stats_close
};

/* ---------------- 约束函数 ---------------- */
/*
 * 这里的全局上下界用于快速排除明显不可能的序列。
 *
 * - global_lower_bound：最小可能的末项位置，确保序列增长不低于“每次至少加 2”的要求。
 * - global_upper_bound_original：从二次方程解出当前前缀时可允许的最大末项值。
 *
 * 若 m 位于 [lb, ub] 之外，则当前状态不可能再构成合法前缀，直接剪枝。
 */
ll global_lower_bound(ll a1, ll n) {
    return a1 + 2 * (n - 1);
}

ll global_upper_bound_original(ll a1, ll Tn) {
    long double a1_minus1 = (long double)(a1 - 1);
    long double disc = a1_minus1 * a1_minus1 + 4.0L * (long double)Tn;
    long double sqrt_disc = sqrtl(disc);
    long double result = (-a1_minus1 + sqrt_disc) / 2.0L;

    ll ub = (ll)floorl(result);
    while ((long double)(ub + 1) * (long double)(ub + 1) +
               a1_minus1 * (long double)(ub + 1) - (long double)Tn <= 0) {
        ub++;
    }
    while ((long double)ub * (long double)ub +
               a1_minus1 * (long double)ub - (long double)Tn > 0) {
        ub--;
    }
    return ub;
}

bool check_refined_ineq(State *s, ll L) {
    ll A = s->m;
    ll n = s->n;
    ll B = s->a[A];
    ll a1 = s->a[1];

    long double left = (long double)a1 +
                       (long double)(n - 1) * (long double)(A - n + 2) +
                       (long double)(A - n) * (long double)(B - A + n + 1);
    long double right = (long double)L * ((long double)B + (long double)L + 1.0L);

    return left >= right;
}

bool check_L_equals_upper_extra(State *s, ll L) {
    if (L != s->m - s->n) return true;
    long double diff = (long double)(s->m - s->n);
    long double threshold = 1.0L + sqrtl((long double)(1 + s->a[1]));
    return diff <= threshold;
}

/* ---------------- 前向声明 ---------------- */
bool fill_segment(State *s, ll old_m, ll L, ll pos, ll need_sum, ll min_val);
bool try_extend(State *s);

/* ---------------- 填充新增段（迭代回溯版） ---------------- */
typedef struct {
    ll pos;
    ll min_val;
    ll max_val;
    ll curr_val;
    ll need_sum;
    ll idx;
} FillFrame;
/*
 * fill_segment：在当前前缀末端插入一个新的区间，长度为 L，
 * 其目标是让新扩展段的元素值满足严格递增且仍保持前缀和约束。
 *
 * 这部分使用“迭代回溯”而不是递归：
 *   - stack 模拟递归调用栈；
 *   - 每个 FillFrame 保存当前扩展位置、最小值、最大值、当前尝试值和剩余和；
 *   - 只要找到一条可行分支，就将新段写回 s->a 并继续尝试扩展。
 *
 * 也就是说，这里是真正的“DFS 但无递归栈爆炸”的实现。
 */bool fill_segment(State *s, ll old_m, ll L, ll pos, ll need_sum, ll min_val) {
    if (L == 0) return (need_sum == 0);

    size_t stack_cap = 16;
    FillFrame *stack = (FillFrame *)safe_realloc(NULL, stack_cap * sizeof(FillFrame), s);

    ll top = 0;
    stack[top].pos = pos;
    stack[top].min_val = min_val;
    stack[top].need_sum = need_sum;
    stack[top].idx = old_m + pos + 1;

    ll prev_idx = stack[top].idx - 1;
    ll prev_val = (pos == 0) ? s->a[old_m] : s->a[prev_idx];
    ll max_allowed = 2 * prev_val - prev_idx;
    if (max_allowed > MAX_VAL) max_allowed = MAX_VAL;

    ll rem = L - pos - 1;
    ll val_max_sum;
    if (rem == 0) {
        val_max_sum = need_sum;
    } else {
        ll tmp;
        if (__builtin_mul_overflow(rem, rem + 1, &tmp)) {
            free(stack);
            return false;
        }
        if (tmp > need_sum) {
            free(stack);
            return false;
        }
        val_max_sum = (need_sum - tmp) / (rem + 1);
    }
    stack[top].max_val = (val_max_sum < max_allowed) ? val_max_sum : max_allowed;
    if (stack[top].max_val < min_val) {
        free(stack);
        return false;
    }
    stack[top].curr_val = min_val - 1;

    while (top >= 0) {
        FillFrame *f = &stack[top];

        if (f->pos == L) {
            if (f->need_sum == 0) {
                ll old_m_saved = s->m;
                ll old_n_saved = s->n;
                ll old_Tn_saved = s->Tn;

                s->m = old_m + L;
                s->n = s->n + 1;
                s->Tn = old_Tn_saved + old_Tn_saved;

                bool ok = try_extend(s);

                s->m = old_m_saved;
                s->n = old_n_saved;
                s->Tn = old_Tn_saved;

                free(stack);
                return ok;
            } else {
                top--;
                continue;
            }
        }

        f->curr_val++;
        if (f->curr_val > f->max_val) {
            top--;
            continue;
        }

        ll val = f->curr_val;
        ensure_capacity(s, f->idx);
        s->a[f->idx] = val;

        ll next_pos = f->pos + 1;
        ll next_need_sum = f->need_sum - val;
        ll next_min_val = val + 2;
        ll next_idx = f->idx + 1;

        if (next_pos == L) {
            if ((size_t)top + 1 >= stack_cap) {
                stack_cap *= 2;
                stack = (FillFrame *)safe_realloc(stack, stack_cap * sizeof(FillFrame), s);
            }
            top++;
            stack[top].pos = L;
            stack[top].need_sum = next_need_sum;
            stack[top].min_val = next_min_val;
            stack[top].max_val = 0;
            stack[top].curr_val = 0;
            stack[top].idx = next_idx;
        } else {
            prev_idx = next_idx - 1;
            prev_val = val;
            max_allowed = 2 * prev_val - prev_idx;
            if (max_allowed > MAX_VAL) max_allowed = MAX_VAL;

            ll next_rem = L - next_pos - 1;
            ll next_val_max_sum;
            if (next_rem == 0) {
                next_val_max_sum = next_need_sum;
            } else {
                ll tmp;
                if (__builtin_mul_overflow(next_rem, next_rem + 1, &tmp)) {
                    continue;
                }
                if (tmp > next_need_sum) {
                    continue;
                }
                next_val_max_sum = (next_need_sum - tmp) / (next_rem + 1);
            }
            ll next_max_allowed_final = (next_val_max_sum < max_allowed) ? next_val_max_sum : max_allowed;
            if (next_max_allowed_final < next_min_val) {
                continue;
            }

            if ((size_t)top + 1 >= stack_cap) {
                stack_cap *= 2;
                stack = (FillFrame *)safe_realloc(stack, stack_cap * sizeof(FillFrame), s);
            }
            top++;
            stack[top].pos = next_pos;
            stack[top].min_val = next_min_val;
            stack[top].max_val = next_max_allowed_final;
            stack[top].curr_val = next_min_val - 1;
            stack[top].need_sum = next_need_sum;
            stack[top].idx = next_idx;
        }
    }

    free(stack);
    return false;
}

/* ---------------- 尝试扩展一步 ---------------- */
/*
 * try_extend 是这段程序的核心入口：
 *   - 先判断当前状态是否已经达到 promising 条件；
 *   - 若不满足，则尝试选取一个新的长度 L，构造下一段；
 *   - 通过几类数学不等式和和的下界过滤，剪去不可能扩展的分支。
 *
 * 这里的策略是：
 *   1. 先做全局范围检查；
 *   2. 再检查已知步长 L 是否合理；
 *   3. 若 L 未确定，则枚举所有可能 L；
 *   4. 对每个 L 调用 fill_segment 进行具体填充。
 */
bool try_extend(State *s) {
    if (s->Tn > LIMIT ||
        s->n > PROMISING_N_THRESHOLD ||
        s->a[s->m] > PROMISING_ITEM_THRESHOLD) {
        log_attempt(s, FAIL_LIMIT_REACHED, "try_extend: 满足 promising 条件");
        save_promising_sequence(s);
        return true;
    }

    ll lb = global_lower_bound(s->a[1], s->n);
    ll ub_orig = global_upper_bound_original(s->a[1], s->Tn);
    if (s->m < lb || s->m > ub_orig) {
        log_attempt(s, FAIL_STEP_ABOVE_MAX, "try_extend: 全局上界/下界检查失败");
        return false;
    }

    ll old_m = s->m;
    ll old_n = s->n;
    ll old_Tn = s->Tn;
    ll L;

    if (old_n + 1 <= old_m) {
        L = s->a[old_n + 1] - s->a[old_n];

        if (L < 2) {
            log_attempt(s, FAIL_STEP_BELOW_MIN, "try_extend: 已知步长 < 2");
            return false;
        }
        if (L > old_m - old_n) {
            log_attempt(s, FAIL_STEP_ABOVE_MAX, "try_extend: 已知步长超过上界");
            return false;
        }
        if (!check_L_equals_upper_extra(s, L)) {
            log_attempt(s, FAIL_L_EXTRA, "try_extend: L取上界额外限制失败");
            return false;
        }
        if (!check_refined_ineq(s, L)) {
            log_attempt(s, FAIL_REFINED_INEQ, "try_extend: 精细不等式不成立");
            return false;
        }

        ll start_min = s->a[old_m] + 2;
        ll part1, part2;
        if (mul_overflow(L, start_min, &part1)) return false;
        if (mul_overflow(L, L - 1, &part2)) return false;
        ll min_sum;
        if (add_overflow(part1, part2, &min_sum)) return false;

        if (min_sum > old_Tn) {
            log_attempt(s, FAIL_MIN_SUM_EXCEED, "try_extend: 新增段最小和超过 Tn");
            return false;
        }

        return fill_segment(s, old_m, L, 0, old_Tn, start_min);

    } else {
        ll upper = old_m - old_n;
        for (L = 2; L <= upper; L++) {
            if (old_m + L > MAX_VAL) break;
            if (!check_refined_ineq(s, L)) {
                log_attempt(s, FAIL_REFINED_INEQ, "try_extend: 枚举L精细不等式失败");
                break;
            }
            if (!check_L_equals_upper_extra(s, L)) {
                log_attempt(s, FAIL_L_EXTRA, "try_extend: 枚举L额外限制失败");
                continue;
            }
            ll start_min = s->a[old_m] + 2;
            ll part1, part2;
            if (mul_overflow(L, start_min, &part1)) break;
            if (mul_overflow(L, L - 1, &part2)) break;
            ll min_sum;
            if (add_overflow(part1, part2, &min_sum)) break;

            if (min_sum > old_Tn) {
                log_attempt(s, FAIL_MIN_SUM_EXCEED, "try_extend: 枚举L最小和超过Tn");
                break;
            }
            if (fill_segment(s, old_m, L, 0, old_Tn, start_min)) {
                return true;
            }
        }
        log_attempt(s, FAIL_FILL_IMPOSSIBLE, "try_extend: 枚举L全部失败");
        return false;
    }
}

/* ---------------- 枚举初始段 ---------------- */
/*
 * enumerate_initial 负责从 a1 开始，构造最初的那一段前缀。
 *
 * 它的逻辑是：
 *   - 固定 a1 = A；
 *   - 枚举 a2, a3, ... 直到需要的长度耗尽；
 *   - 每个新值都必须满足“最大允许值”约束 2 * prev - prev_idx；
 *   - 一旦前缀形成合法状态，就把控制交给 try_extend 继续扩展。
 *
 * 换句话说，这一层负责把所有候选初始前缀穷举出来，
 * 再交由后续判定筛出真正可扩展的 promising 序列。
 */
bool enumerate_initial(State *s, ll A, ll cur_idx, ll need_terms, ll min_val, ll sum_so_far) {
    if (need_terms == 0) {
        s->m = A;
        s->n = 1;
        s->Tn = sum_so_far;

        long double a1_minus1_sq = (long double)(A - 1) * (A - 1);
        long double term = (a1_minus1_sq - 1.0L) * (a1_minus1_sq - 1.0L) /
                           (16.0L * (long double)sum_so_far);
        if (term <= 0.0L) {
            s->N_threshold = 1;
        } else {
            s->N_threshold = (ll)ceill(1.0L + log2l(term));
        }

        return try_extend(s);
    }

    ll prev = s->a[cur_idx - 1];
    ll prev_idx = cur_idx - 1;
    ll max_allowed = 2 * prev - prev_idx;
    if (max_allowed > MAX_VAL) max_allowed = MAX_VAL;
    if (max_allowed < min_val) {
        log_attempt(s, FAIL_INITIAL_RANGE_EMPTY, "enumerate_initial: 枚举范围为空");
        return false;
    }

    for (ll val = min_val; val <= max_allowed; val++) {
        ensure_capacity(s, cur_idx);
        s->a[cur_idx] = val;
        if (enumerate_initial(s, A, cur_idx + 1, need_terms - 1, val + 2, sum_so_far + val)) {
            return true;
        }
    }

    log_attempt(s, FAIL_INITIAL_RANGE_EMPTY, "enumerate_initial: 所有val失败");
    return false;
}

/* ---------------- 主函数 ---------------- */
/*
 * 主循环的查找过程：
 *   1. 固定 a1 的值，从 2 到 max_a1 逐个枚举；
 *   2. 对每个 a1，枚举所有可能的 a2；
 *   3. 再调用 enumerate_initial 构造初始前缀；
 *   4. 只要找到一个可扩展的 promising 序列，就输出并结束。
 *
 * 这对应“按起点枚举 -> 按二项构造 -> 逐步搜索扩展”的整体思路。
 */
int main(void) {
    // 默认使用详细文件日志
    if (!init_detailed_log("r2.log")) {
        fprintf(stderr, "错误：无法打开日志文件 r2.log\n");
        return 1;
    }
    // 切换日志模式：
    // current_sink = &stats_sink;   // 统计日志
    // current_sink = &null_sink;    // 无日志

    ll max_a1 = 50;

    for (ll A = 2; A <= max_a1; A++) {
        printf("正在测试 a1 = %lld ...\n", A);

        st.cap = 0;
        st.a = NULL;
        ensure_capacity(&st, A);
        st.a[1] = A;

        ll a2_min = A + 2;
        ll a2_max = 2 * A - 2;
        if (a2_max > MAX_VAL) a2_max = MAX_VAL;

        bool found = false;
        for (ll a2 = a2_min; a2 <= a2_max; a2++) {
            st.a[2] = a2;
            ll sum = A + a2;

            if (enumerate_initial(&st, A, 3, A - 2, a2 + 2, sum)) {
                found = true;
                break;
            }
        }

        if (found) {
            printf("发现一个promising序列: a1=%lld\n", A);
            if (promising_a && promising_m > 0) {
                printf("PROMISING: a_1..a_%lld =", promising_m);
                ll limit = promising_m < 200 ? promising_m : 200;
                for (ll i = 1; i <= limit; i++) {
                    printf(" %lld", promising_a[i]);
                }
                if (promising_m > 200) printf(" ... (共 %lld 项)", promising_m);
                printf("\n");
            }
        }

        free(st.a);
        st.a = NULL;

        if (found) {
            break;
        }
    }

    current_sink->flush();
    current_sink->close();

    if (promising_a) free(promising_a);

    return 0;
}
