# Todo list

[Frontend](/oss/python/deepagents/frontend/overview)[Patterns](/oss/python/deepagents/frontend/subagent-streaming)

# Todo list
Copy page

Track agent progress with a real-time todo list synced from agent stateCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Not every agent interaction is a chat. Sometimes the agent is executing a
multi-step plan, and the best way to show progress is a **todo list** that
updates in real time. The deep agent todo list pattern reads a `todos` array
directly from the agent’s state, rendering each item with its current status as
the agent works through its plan. It’s a progress dashboard built on the same
`useStream` hook you use for chat. It shows that agent state can power any UI,
not just message bubbles.

## [​](#how-it-works)How it works

In a LangGraph agent, state isn’t limited to messages. You can define **custom
state keys** that hold arbitrary data. In this case, a `todos` array. As the
agent executes its plan, it updates each todo’s status from `"pending"` to
`"in_progress"` to `"completed"`. The `useStream` hook exposes these custom
state values via `stream.values`, and your UI renders them reactively.
The flow looks like this:

{' ' * (self.list_depth - 1)}- User submits a request

{' ' * (self.list_depth - 1)}- Agent creates a plan and populates `todos` in its state

{' ' * (self.list_depth - 1)}- Agent begins executing each todo transitions through `pending` →
`in_progress` → `completed`

{' ' * (self.list_depth - 1)}- `stream.values.todos` updates in real time as the agent progresses

{' ' * (self.list_depth - 1)}- Your UI re-renders the todo list with current statuses

## [​](#setting-up-usestream)Setting up useStream

No special configuration is needed. Point `useStream` at your agent and
read the `todos` from `stream.values`.
Define a TypeScript interface matching your agent’s state schema and pass it as a type parameter to `useStream` for type-safe access to state values, including custom state keys like `todos`. In the examples below, replace `typeof myAgent` with your interface name:

```
import type { BaseMessage } from "@langchain/core/messages";

interface TodoItem {
 title: string;
 status: "pending" | "in_progress" | "completed";
 description?: string;
}

interface AgentState {
 messages: BaseMessage[];
 todos: TodoItem[];
}

```

ReactVueSvelteAngular
```
import { useStream } from "@langchain/react";

const AGENT_URL = "http://localhost:2024";

export function TodoAgent() {
 const stream = useStream<typeof myAgent>({
 apiUrl: AGENT_URL,
 assistantId: "deep_agent_todo_list",
 });

 const todos = stream.values?.todos ?? [];

 return (
 <div>
 <TodoList todos={todos} />
 {stream.messages.map((msg) => (
 <Message key={msg.id} message={msg} />
 ))}
 </div>
 );
}

```

## [​](#the-todo-interface)The Todo interface

Each todo in the array has a simple structure:

```
interface Todo {
 status: "pending" | "in_progress" | "completed";
 content: string;
}

```

| 
 | Property | Description
 | `status` | The current state of this task. Options: `pending` (not started), `in_progress` (agent is working on it), `completed` (done)
 | `content` | A human-readable description of what the task involves
The agent populates this array when it creates its plan, then updates individual
items as it executes each step.

## [​](#building-the-todolist-component)Building the TodoList component

The todo list renders each item with a status icon, color coding, and visual
styling that reflects the current state:

```
function TodoList({ todos }: { todos: Todo[] }) {
 const completed = todos.filter((t) => t.status === "completed").length;
 const percentage = todos.length
 ? Math.round((completed / todos.length) * 100)
 : 0;

 return (
 <div className="rounded-lg border bg-white p-4 shadow-sm">
 <div className="mb-4 flex items-center justify-between">
 <h2 className="text-lg font-semibold">Agent Progress</h2>
 <span className="text-sm text-gray-500">
 {completed}/{todos.length} tasks
 </span>
 </div>

 <ProgressBar percentage={percentage} />

 <ul className="mt-4 space-y-2">
 {todos.map((todo, i) => (
 <TodoItem key={i} todo={todo} />
 ))}
 </ul>
 </div>
 );
}

```

## [​](#progress-bar)Progress bar

A visual progress bar gives users an at-a-glance summary of overall completion:

```
function ProgressBar({ percentage }: { percentage: number }) {
 return (
 <div className="space-y-1">
 <div className="flex items-center justify-between text-xs text-gray-500">
 <span>Progress</span>
 <span>{percentage}%</span>
 </div>
 <div className="h-2 overflow-hidden rounded-full bg-gray-200">
 <div
 className="h-full rounded-full bg-green-500 transition-all duration-500"
 style={{ width: `${percentage}%` }}
 />
 </div>
 </div>
 );
}

```

## [​](#individual-todo-items)Individual todo items

Each item gets a status icon, color-coded text, and strikethrough styling for
completed tasks:

```
function TodoItem({ todo }: { todo: Todo }) {
 const config = {
 pending: {
 icon: "○",
 textClass: "text-gray-600",
 bgClass: "bg-gray-50",
 iconClass: "text-gray-400",
 },
 in_progress: {
 icon: "◉",
 textClass: "text-amber-800",
 bgClass: "bg-amber-50 border-amber-200",
 iconClass: "text-amber-500 animate-pulse",
 },
 completed: {
 icon: "✓",
 textClass: "text-green-800 line-through",
 bgClass: "bg-green-50 border-green-200",
 iconClass: "text-green-500",
 },
 };

 const style = config[todo.status];

 return (
 <li
 className={`flex items-start gap-3 rounded-md border px-3 py-2 ${style.bgClass}`}
 >
 <span className={`mt-0.5 text-lg leading-none ${style.iconClass}`}>
 {style.icon}
 </span>
 <span className={`text-sm ${style.textClass}`}>{todo.content}</span>
 </li>
 );
}

```

The `in_progress` icon uses `animate-pulse` to draw attention to the currently
active task.

## [​](#calculating-progress)Calculating progress

Derive progress metrics directly from the todos array:

```
const todos = stream.values?.todos ?? [];

const completed = todos.filter((t) => t.status === "completed").length;
const inProgress = todos.filter((t) => t.status === "in_progress").length;
const pending = todos.filter((t) => t.status === "pending").length;
const percentage = todos.length
 ? Math.round((completed / todos.length) * 100)
 : 0;

```

These values update reactively as the agent modifies its state, keeping the
progress bar and counters in sync.

## [​](#combining-with-chat-messages)Combining with chat messages

The todo list works alongside the regular chat interface. A practical layout
shows the todo list as a persistent sidebar or header panel, with chat messages
below:

```
function TodoAgentLayout() {
 const stream = useStream<typeof myAgent>({
 apiUrl: AGENT_URL,
 assistantId: "deep_agent_todo_list",
 });

 const todos = stream.values?.todos ?? [];

 return (
 <div className="flex h-screen flex-col">
 {todos.length > 0 && (
 <div className="border-b bg-gray-50 p-4">
 <TodoList todos={todos} />
 </div>
 )}

 <main className="flex-1 overflow-y-auto p-6">
 <div className="mx-auto max-w-2xl space-y-4">
 {stream.messages.map((msg) => (
 <Message key={msg.id} message={msg} />
 ))}
 </div>
 </main>

 <ChatInput
 onSubmit={(text) =>
 stream.submit({ messages: [{ type: "human", content: text }] })
 }
 isLoading={stream.isLoading}
 />
 </div>
 );
}

```

Show the todo list only when `todos.length > 0`. Before the agent creates its
plan, there’s nothing to display. Showing an empty component wastes space.

## [​](#custom-state-beyond-todos)Custom state beyond todos

This pattern demonstrates a powerful principle: `stream.values` can expose
**any custom state** your agent defines, not just messages. The `todos` array is
just one example. You could use the same approach for:

{' ' * (self.list_depth - 1)}- **Progress metrics**: `stream.values.progress` with numeric completion data

{' ' * (self.list_depth - 1)}- **Generated artifacts**: `stream.values.document` with a structured document
the agent is building

{' ' * (self.list_depth - 1)}- **Decision logs**: `stream.values.decisions` tracking every choice the agent
made

{' ' * (self.list_depth - 1)}- **Resource lists**: `stream.values.sources` with links and references the
agent found

```
// Any custom state key your agent defines is accessible
const document = stream.values?.document;
const sources = stream.values?.sources ?? [];
const confidence = stream.values?.confidence_score;

```

Custom state keys are defined in your LangGraph graph’s state schema. The
`useStream` hook automatically includes them in `stream.values` without any additional
client-side configuration.

## [​](#animating-transitions)Animating transitions

Todo status transitions happen in real time, and smooth animations make these
changes feel polished rather than jarring:

```
function TodoItem({ todo }: { todo: Todo }) {
 return (
 <li
 className={`
 flex items-start gap-3 rounded-md border px-3 py-2
 transition-all duration-300 ease-in-out
 ${getStatusStyles(todo.status)}
 `}
 >
 <span
 className={`
 mt-0.5 text-lg leading-none transition-colors duration-300
 ${getIconStyles(todo.status)}
 `}
 >
 {getStatusIcon(todo.status)}
 </span>
 <span
 className={`
 text-sm transition-all duration-300
 ${todo.status === "completed" ? "line-through opacity-60" : ""}
 `}
 >
 {todo.content}
 </span>
 </li>
 );
}

```

The `transition-all duration-300` classes ensure that color changes,
strikethrough, and opacity shifts all animate smoothly.

## [​](#use-cases)Use cases

The todo list pattern fits any scenario where an agent executes a structured
plan:

{' ' * (self.list_depth - 1)}- **Project planning**: agent breaks a project into tasks and works through
them sequentially

{' ' * (self.list_depth - 1)}- **Research workflows**: each research question becomes a todo that the agent
investigates and completes

{' ' * (self.list_depth - 1)}- **Data processing**: steps like ingestion, validation, transformation, and
export each get their own todo

{' ' * (self.list_depth - 1)}- **Onboarding flows**: agent walks through setup steps, checking off each one
as it configures services

{' ' * (self.list_depth - 1)}- **Report generation**: sections of a report become todos: gather data,
analyze trends, write summary, format output

## [​](#handling-empty-and-loading-states)Handling empty and loading states

Handle the initial state before the agent has created its plan:

```
function TodoList({ todos, isLoading }: { todos: Todo[]; isLoading: boolean }) {
 if (todos.length === 0 && !isLoading) {
 return null;
 }

 if (todos.length === 0 && isLoading) {
 return (
 <div className="rounded-lg border bg-white p-4 shadow-sm">
 <div className="flex items-center gap-2 text-sm text-gray-500">
 <span className="animate-spin">⟳</span>
 Agent is creating a plan...
 </div>
 </div>
 );
 }

 return (
 <div className="rounded-lg border bg-white p-4 shadow-sm">
 {/* ... full todo list rendering */}
 </div>
 );
}

```

## [​](#best-practices)Best practices

{' ' * (self.list_depth - 1)}- **Show the todo list prominently**. It’s the primary progress indicator for
plan-based agents. Don’t bury it below the fold.

{' ' * (self.list_depth - 1)}- **Animate status transitions**. Smooth transitions make the agent feel more
responsive. Use CSS transitions on background color, text decoration, and
opacity.

{' ' * (self.list_depth - 1)}- **Only highlight one `in_progress` item**. Agents typically work on one task
at a time. If multiple items show as `in_progress`, the UI gets noisy.
Consider only pulsing the first one.

{' ' * (self.list_depth - 1)}- **Collapse or dim completed items**. As the list grows, completed items
become less relevant. Reduce their visual weight so users focus on what’s
still happening.

{' ' * (self.list_depth - 1)}- **Show the progress percentage**. A single number like “67% complete” is
immediately understandable, even from across the room.

{' ' * (self.list_depth - 1)}- **Keep the todo list in sync**. Because `stream.values` updates reactively,
the todo list stays current automatically. Don’t add manual polling or
refresh logic.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/deepagents/frontend/todo-list.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Subagent streamingPrevious](/oss/python/deepagents/frontend/subagent-streaming)[SandboxNext](/oss/python/deepagents/frontend/sandbox)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
